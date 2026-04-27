// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// aflow_pocc.h and aflow_pocc*.cpp*
//
// completely revised approach to KESONG YANG's original implementation
// no recursion needed --- too memory intensive
// issues with UFF (structure comparison), multiply occupied sites, and
// vacancies are all addressed robustly here
// 2013-2019: corey.oses@duke.edu
// 2010-2011: kesong.yang@gmail.com (LEGACY)

#ifndef _AFLOW_POCC_H_
#define _AFLOW_POCC_H_

#include <cassert>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <flow/aflow_pflow.h>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_type_traits.tpp"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xoption.h"
#include "AUROSTD/aurostd_xparser_json.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflow.h"
#include "flow/aflow_xclasses.h"
#include "modules/APL/aflow_apl.h"
#include "structure/aflow_xatom.h"
#include "structure/aflow_xstructure.h"

// precision defines tol and paddings
// const int _AFLOW_POCC_PRECISION_ = 8;    //moved to aflow.h
// tolerances
// const double _AFLOW_POCC_ZERO_TOL_ = pow(10,-_AFLOW_POCC_PRECISION_);  //moved to aflow.h

// uff param modes
const uint BOND_MODE = 0;
const uint NONBOND_MODE = 1;

// files
const std::string AFLOW_POCC_TAG = "[AFLOW_POCC]";

// elements
// to make POCC chemistry-independent, calculate UFF based on standard set of elements
// we should not need to have more than 30 species, but if you want to add to list
// I recommend adding IN ORDER so that elements become more different left to right (like periodic table)
const std::string STD_ELEMENTS_LIST =
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
    "La Hf Ta W Re Os Ir Pt Au Hg ";

const int TEMPERATURE_PRECISION = 2;  // not really going to explore more than 2000-3000K, looks weird if decimal is larger than non-decimal part of number //4;  //this is std::fixed

namespace pocc {
  class POccCalculatorBuilder;
  class POccCalculator; // forward declaration
  class POccSuperCellSet; // forward declaration

  std::string POCC_MINIMUM_CONFIGURATION(const aurostd::xoption& vpflow);
  std::string POCC_MINIMUM_CONFIGURATION(const std::string& directory = "./");
  bool structuresGenerated(const std::string& directory = ".");
  xstructure extractPARTCAR(const std::string& AflowIn);

  std::vector<double> getVTemperatures(const std::string& temp_string);

  void parsePOccHashFromXStructureTitle(const std::string& title, std::string& pocc_hash);
  void parsePOccHashFromXStructureTitle(const std::string& title, std::string& pocc_hash, std::string& hnf_index_str, std::string& site_config_index_str);
  unsigned long long int getDGFromXStructureTitle(const std::string& title);
  void parsePropertyByTag(const std::string& line, const std::string& tag, double& prop);
  bool patchStructuresAllFile(const _aflags& aflags, std::string& structures_file, std::stringstream& structures_file_ss, std::ofstream& FileMESSAGE, std::ostream& oss);
  bool patchStructuresAllFile(std::stringstream& structures_file_ss);
  bool patchStructuresUniqueFile(const _aflags& aflags, std::string& structures_file, std::stringstream& structures_file_ss, std::ofstream& FileMESSAGE, std::ostream& oss);
  bool patchStructuresUniqueFile(std::stringstream& structures_file_ss);
  void patchStructuresFile(const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);

  void updateProgressBar(unsigned long long int current, unsigned long long int end, std::ostream& oss = std::cout);
  std::vector<std::string> getElementsList();

  class POccUnit;  // forward declaration
  class POccSiteConfiguration;    // forward declaration
  class UFFParamAtom;  // forward declaration
  struct UFFParamBond;  // forward declaration

  // bool sortAtoms(const _atom& a1,const _atom& a2);
  // bool getNextSiteConfiguration(vector<vector<POccSiteConfiguration> >& vv_count_configs,vector<uint>& v_config_order,vector<int>& v_config_iterators,vector<vector<int> >& v_types_config);
  bool getNextSiteConfiguration(std::vector<std::vector<POccSiteConfiguration>>& vv_count_configs, std::vector<int>& v_config_iterators, std::vector<std::vector<int>>& v_types_config);
  bool getNextSiteConfiguration(std::vector<uint>& v_config_order, std::vector<std::vector<int>>& v_types_config);
  bool getNextSiteConfiguration(std::vector<int>& site_config);
  std::vector<uint> getConfigOrder(std::vector<std::vector<int>>& v_types_config);
  std::string getUFFParameterString(const std::string& element);
  std::vector<UFFParamAtom> getTypes2UFFParamsMap(const std::vector<std::string>& elements);
  xstructure createNonPOccStructure(const xstructure& xstr_pocc, const std::vector<POccUnit>& pocc_sites, bool use_automatic_volumes_in = true);             // convert pocc xstructure to non-pocc
  std::vector<POccUnit> getPOccSites(const xstructure& xstr_pocc, std::ostream& oss = std::cout);
  std::vector<POccUnit> getPOccSites(const xstructure& xstr_pocc, const _aflags& aflags, std::ostream& oss = std::cout);
  std::vector<POccUnit> getPOccSites(const xstructure& xstr_pocc, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  std::vector<POccUnit> getPOccSites(const xstructure& xstr_pocc, const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
  std::vector<uint> getTypes2PCMap(const xstructure& xstr);

  bool sortPOccSites(const POccUnit& p1, const POccUnit& p2);
  bool sortPOccGroups(const POccUnit& p1, const POccUnit& p2);

  std::string getARUNString(const std::list<POccSuperCellSet>& l_supercell_sets, unsigned long long int i);
  std::string getARUNString(unsigned long long int index_structure_group,
                            unsigned long long int vstructure_groups_size,
                            unsigned long long int index_structure,
                            unsigned long long int vstructures_size,
                            unsigned long long int index_hnf,
                            unsigned long long int index_site_config,
                            bool include_strgrp = false);

  double getHmix(const aurostd::xvector<double>& xv_energies, const aurostd::xvector<double>& xv_dgs);
  double getHmix(const aurostd::xvector<double>& xv_energies, const aurostd::xvector<double>& xv_dgs, double& dg_total);
  double getEFA(const aurostd::xvector<double>& xv_energies, const aurostd::xvector<double>& xv_dgs);

  void poccOld2New(std::ostream& oss = std::cout);
  void poccOld2New(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);

  void POCC_Convolution(const aurostd::xoption& vpflow);

  std::string addDefaultPOCCTOL2string(const std::string& input);

  void getTemperatureStringParameters(int& temperature_precision, bool& temperatures_int, int& zero_padding_temperature);
  void getTemperatureStringParameters(const std::vector<double>& v_temperatures, int& temperature_precision, bool& temperatures_int, int& zero_padding_temperature);
  std::string getTemperatureString(double temperature, int precision, bool temperatures_int, int zero_padding_temperature);

  double poccDOSCAR2temperature(const std::string& doscar_path);
} // namespace pocc

namespace pocc {
  // POccUnit could be
  // 1) a bunch of atoms on the same site (site)
  // 2) within a site, a bunch of atoms with the same partial_occupation_value (group)
  class POccUnit : public xStream, public JsonSerializable<POccUnit> {
  public:
    POccUnit(std::ostream& oss = std::cout) : xStream(oss) {};
    POccUnit(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss) {}
    POccUnit(const _aflags& aflags, std::ostream& oss = std::cout) : xStream(oss), m_aflags(aflags) {}
    POccUnit(const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss), m_aflags(aflags) {}
      // bool operator<(const POccUnit& other) const;
    void clear();

    bool m_initialized = false;
    _aflags m_aflags;
    uint site = AUROSTD_MAX_UINT; // reflexive
    bool partial_occupation_flag = false;
    double partial_occupation_value = AUROSTD_MAX_DOUBLE;
    std::vector<uint> v_occupants;
    std::vector<uint> v_types;
    std::vector<POccUnit> m_pocc_groups;
    bool is_inequivalent = false; // from iatoms of non-pocc structure
    uint equivalent = AUROSTD_MAX_UINT; // from iatoms of non-pocc structure

    void setAFlags(const _aflags& aflags);

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccUnit deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccUnit"; }

    // SERIALIZATION MEMBERS
#define JSON_POccUnit_MEMBERS m_initialized, m_aflags, site, partial_occupation_flag, partial_occupation_value, v_occupants, v_types, m_pocc_groups, is_inequivalent, equivalent
  };
}  // namespace pocc

namespace pocc {
  // simple structure for sorting by vacancy count
  struct SiteVacancyCount {
    uint site;
    uint vacancy_count;
    bool operator<(const SiteVacancyCount& other) const;
  };

  struct StructureConfiguration : public JsonSerializable<StructureConfiguration> {
    std::vector<POccSiteConfiguration> site_configs;
    double max_stoich_error;
    double max_site_error;

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    StructureConfiguration deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "StructureConfiguration"; }

    // SERIALIZATION MEMBERS
#define JSON_StructureConfiguration_MEMBERS site_configs, max_stoich_error, max_site_error
  };
} // namespace pocc

namespace pocc {
  class POccSuperCell : public JsonSerializable<POccSuperCell> {
  public:
    POccSuperCell() = default;
    bool operator<(const POccSuperCell& other) const;
    void clear();

    unsigned long long int m_hnf_index = AUROSTD_MAX_ULLINT;
    unsigned long long int m_site_config_index = AUROSTD_MAX_ULLINT;
    unsigned long long int m_degeneracy = 0; // VERY IMPORTANT, 0 m_degeneracy will kill division
    double m_energy_uff = AUROSTD_MAX_DOUBLE;

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccSuperCell deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccSuperCell"; }

    // SERIALIZATION MEMBERS
#define JSON_POccSuperCell_MEMBERS m_hnf_index, m_site_config_index, m_degeneracy, m_energy_uff
  };
  bool sortPSCsUFFEnergy(const POccSuperCell& a, const POccSuperCell& b);
} // namespace pocc

namespace pocc {
  class POccSuperCellSet : public JsonSerializable<POccSuperCellSet> {
  public:
    POccSuperCellSet() = default;
    bool operator<(const POccSuperCellSet& other) const;
    void clear();

    std::vector<POccSuperCell> m_psc_set;
    double m_energy_dft = AUROSTD_MAX_DOUBLE;  // only calculate this for the set
    double m_probability = AUROSTD_MAX_DOUBLE;

    [[nodiscard]] unsigned long long int getDegeneracy() const;
    [[nodiscard]] const POccSuperCell& getSuperCell() const;
    [[nodiscard]] double getHNFIndex() const; // ME20211006
    [[nodiscard]] aurostd::xmatrix<int> getHNFMatrix() const;
    [[nodiscard]] double getSiteConfigIndex() const;
    [[nodiscard]] std::vector<std::vector<int>> getSiteConfig() const; // ME20211006
    [[nodiscard]] xstructure getStructure() const; // ME20211006
    [[nodiscard]] double getUFFEnergy() const;

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccSuperCellSet deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccSuperCellSet"; }

    // SERIALIZATION MEMBERS
#define JSON_POccSuperCellSet_MEMBERS m_psc_set, m_energy_dft, m_probability
  };
  bool sortPSCSetsUFFEnergy(const POccSuperCellSet& a, const POccSuperCellSet& b);
} // namespace pocc

namespace pocc {
  class UFFParamAtom : public JsonSerializable<UFFParamAtom> {
  public:
    UFFParamAtom() = default;
    void clear() { *this = UFFParamAtom(); }

    std::string symbol;
    double r1 = AUROSTD_MAX_DOUBLE;        // bond distance
    double theta0 = AUROSTD_MAX_DOUBLE;
    double x1 = AUROSTD_MAX_DOUBLE;        // nonbond distance
    double D1 = AUROSTD_MAX_DOUBLE;        // nonbond energy
    double zeta = AUROSTD_MAX_DOUBLE;      // scale
    double Z1 = AUROSTD_MAX_DOUBLE;        // effective charge
    double Vi = AUROSTD_MAX_DOUBLE;
    double Uj = AUROSTD_MAX_DOUBLE;
    double ChiI = AUROSTD_MAX_DOUBLE;       // electronegativity
    double hard = AUROSTD_MAX_DOUBLE;
    double radius = AUROSTD_MAX_DOUBLE;

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    UFFParamAtom deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "UFFParamAtom"; }

    // SERIALIZATION MEMBERS
#define JSON_UFFParamAtom_MEMBERS symbol, r1, theta0, x1, D1, zeta, Z1, Vi, Uj, ChiI, hard, radius
  };

  class UFFParamBond : public JsonSerializable<UFFParamBond> {
  public:
    void calculate(UFFParamAtom& uffp1, UFFParamAtom& uffp2, double distij);

    double ren;
    double R0;
    double Kij;
    double Xij;
    double Dij;
    double delta;
    double X6;
    double X12;

  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    UFFParamBond deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "UFFParamBond"; }

// SERIALIZATION MEMBERS
#define JSON_UFFParamBond_MEMBERS ren, R0, Kij, Xij, Dij, delta, X6, X12
  };
} // namespace pocc

namespace pocc {
  class POccSiteConfiguration : public JsonSerializable<POccSiteConfiguration> {
  public:
    POccSiteConfiguration() = default;
    POccSiteConfiguration(int _site, int _i_hnf, std::vector<POccUnit>& pocc_groups) : site(_site), i_hnf(_i_hnf), m_pocc_groups(pocc_groups) {}
    void clear();

    int site = 0;                                         // reflexive
    int i_hnf = 0;                                        // reflexive
    bool partial_occupation_flag = false;
      // any vector or aurostd::xvector is over pocc_groups
    std::vector<POccUnit> m_pocc_groups;                  // reflexive
    aurostd::xvector<int> xv_occupation_count_input;      // pre multiplication, from xstr_pocc
    aurostd::xvector<int> xv_occupation_multiple;              // increment with each pocc_group
    aurostd::xvector<int> xv_occupation_count_supercell;           // post multiplication with multiple
    aurostd::xvector<double> xv_partial_occupation_value;
    aurostd::xvector<double> xv_site_error;
      // sum of occupation_count_total and vacancy_count yields i_hnf (total sites)
    int occupation_count_total = 0;
    int vacancy_count = 0;
    double max_site_error = 0.0;

    void prepareNoPOccConfig();
    void preparePOccConfig();
    int getNextOccupationMultiple(int i_hnf, aurostd::xvector<int>& xv_occupation_multiple);
    int calculateOccupationCountTotal(aurostd::xvector<int>& xv_next_occupation_multiple) const;
    void updateOccupationCounts(int _i_hnf, aurostd::xvector<int>& xv_next_occupation_multiple);
    void calculateError();
    [[nodiscard]] bool isPartiallyOccupied() const;
    [[nodiscard]] std::vector<int> getStartingTypesConfiguration() const;

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccSiteConfiguration deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccSiteConfiguration"; }

    // SERIALIZATION MEMBERS
#define JSON_POccSiteConfiguration_MEMBERS                                                                                                                                                                   \
  site, i_hnf, partial_occupation_flag, m_pocc_groups, xv_occupation_count_input, xv_occupation_multiple, xv_occupation_count_supercell, xv_partial_occupation_value, xv_site_error, occupation_count_total, \
      vacancy_count, max_site_error
  };
} // namespace pocc

namespace pocc {
  class POccCalculatorTemplate {
  public:
    POccCalculatorTemplate() = default;

    xstructure xstr_pocc;                   // input from PARTCAR
    aurostd::xoption m_p_flags;             // e.g., vpflow
    _aflags m_aflags;                       // standard aflow flags
    aurostd::xvector<double> stoich_each_type;       // converting deque<double> to aurostd::xvector<double>
    xstructure xstr_nopocc;                 // will contain symmetry objects (_sym_op, pgroups most important here)
    std::vector<uint> types2pc_map;              // list of atom indices where types2pc_map(0) is 1st type 0 atom, types2pc_map(1) is 1st type 1 atom, etc.
    std::vector<std::string> m_species_redecoration;  // species used to redecorate xstr's for consistent symmetry/uff energies

    virtual void setPOccFlags(const aurostd::xoption& pocc_flags);
    virtual void setAFlags(const _aflags& Aflags);                      // standard _aflags
    virtual void setPOccStructure(const xstructure& xstr_pocc);
    void setNonPOccStructure(const xstructure& xstr_nonpocc);
    virtual void setSpeciesRedecoration(const std::vector<std::string>& species_redecoration);

    // SERIALIZATION MEMBERS, not serializing xoption vars
#define JSON_POccCalculatorTemplate_MEMBERS xstr_pocc, m_aflags, stoich_each_type, xstr_nopocc, types2pc_map, m_species_redecoration
  };
} // namespace pocc

namespace pocc {
  std::vector<uint> getVacanciesSuperCell(const std::vector<int>& pc2sc_map, const std::vector<std::vector<int>>& v_types_config);
  void replaceRandomSitesSuperCell(const xstructure& xstr_pocc, const std::vector<uint>& types2pc_map, const std::vector<int>& pc2sc_map, const std::vector<std::vector<int>>& v_types_config, xstructure& supercell);
  void rebuildSuperCell(const xstructure& xstr_pocc, const std::vector<uint>& v_vacancies, xstructure& supercell);
} // namespace pocc

namespace pocc {
  class POccUFFEnergyAnalyzer;
  /// @brief Class to ease the construction of @c POccUFFEnergyAnalyzer instances.
  /// Simple builder pattern for currying.
  /// Use the @c with method to chain build arguments and then call @c build to get the constructed @c POccUFFEnergyAnalyzer instance.
  /// For example:
  /// @code POccUFFEnergyAnalyzer analyzer = POccUFFEnergyAnalyzerBuilder().with(aflags).with(pocc_flags).build(); @endcode
  ///
  /// @authors
  /// @xlink{https://en.wikipedia.org/wiki/Builder_pattern}
  /// @mod{ST,20250717,created}
  class POccUFFEnergyAnalyzerBuilder {
  public:
    explicit POccUFFEnergyAnalyzerBuilder(std::ostream& oss = std::cout) : p_oss(&oss) {}
    explicit POccUFFEnergyAnalyzerBuilder(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : p_oss(&oss), p_ofs(&FileMESSAGE) {}

    POccUFFEnergyAnalyzerBuilder& with(const std::pair<const xstructure&, const xstructure&>& xstrs) {
      p_xstr_pocc = &xstrs.first;
      p_xstr_nopocc = &xstrs.second;
      return *this;
    }
    POccUFFEnergyAnalyzerBuilder& with(const aurostd::xoption& pocc_flags) {
      p_pocc_flags = &pocc_flags;
      return *this;
    }
    POccUFFEnergyAnalyzerBuilder& with(const _aflags& aflags) {
      p_aflags = &aflags;
      return *this;
    }
    POccUFFEnergyAnalyzerBuilder& with(const std::vector<std::string>& species_redecoration) {
      p_species_redecoration = &species_redecoration;
      return *this;
    }
    POccUFFEnergyAnalyzerBuilder& with(std::ostream& oss) {
      p_oss = &oss;
      return *this;
    }
    POccUFFEnergyAnalyzerBuilder& with(std::ofstream& ofs) {
      p_ofs = &ofs;
      return *this;
    }
    template <typename... BuildArgs> POccUFFEnergyAnalyzerBuilder& with_all(BuildArgs&&... args) {
      (with(std::forward<BuildArgs>(args)), ...);
      return *this;
    }
    [[nodiscard]] POccUFFEnergyAnalyzer build() const;

    // ST20250715 note: the `with` overloads only work if the arguments are all different types
    // the two xstructures are handled by making them a pair, but they need to be given as a pair

  private:
    std::ostream* p_oss = nullptr;
    std::ofstream* p_ofs = nullptr;
    const _aflags* p_aflags = nullptr;
    const xstructure* p_xstr_pocc = nullptr;
    const xstructure* p_xstr_nopocc = nullptr;
    const aurostd::xoption* p_pocc_flags = nullptr;
    const std::vector<std::string>* p_species_redecoration = nullptr;
  };

  class POccUFFEnergyAnalyzer : public POccCalculatorTemplate, public xStream, public JsonSerializable<POccUFFEnergyAnalyzer> {
    friend POccUFFEnergyAnalyzerBuilder;

  public:
    explicit POccUFFEnergyAnalyzer(std::ostream& oss = std::cout) : xStream(oss) {}
    explicit POccUFFEnergyAnalyzer(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss) {}

    template <typename... BuildArgs, aurostd::enable_if_none_of<POccUFFEnergyAnalyzer, BuildArgs...> = true>
    explicit POccUFFEnergyAnalyzer(BuildArgs&&... args) : POccUFFEnergyAnalyzer(build_from(std::forward<BuildArgs>(args)...)) {}
    template <typename... BuildArgs> static POccUFFEnergyAnalyzer build_from(BuildArgs&&... args) { return POccUFFEnergyAnalyzerBuilder().with_all(std::forward<BuildArgs>(args)...).build(); }

    void clear();

      // underlying data structures
    bool m_initialized = false;
    aurostd::xmatrix<double> hnf_mat;
    std::vector<std::vector<int>> m_types_config;            // the config for which we determined bonding
    std::vector<UFFParamAtom> types2uffparams_map;
    std::vector<uint> m_vacancies;
    double m_exploration_radius = AUROSTD_MAX_DOUBLE;
    aurostd::xmatrix<double> distance_matrix;                    // references xstr_nopocc
    std::vector<double> v_dist_nn;                           // references xstr_nopocc
    xstructure xstr_ss;       // superstructure
    std::vector<int> sc2pc_map;
    std::vector<int> pc2sc_map;

    void setSpeciesRedecoration(const std::vector<std::string>& species_redecoration) override;
    void setExplorationRadius();
    void getCluster(aurostd::xmatrix<double>& _hnf_mat);
    void setBonds(const std::vector<std::vector<int>>& v_types_config);
    double getUFFEnergy();
    bool isVacancy(std::vector<uint>& v_vacancies, uint atom);
    double bondEnergyBond(const UFFParamBond& uffb);
    double bondEnergyNoBond(const UFFParamBond& uffb);
    double getUFFBondEnergy(xstructure& xstr, std::vector<std::vector<uint>>& v_bonded_atom_indices, uint MODE);

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccUFFEnergyAnalyzer deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccUFFEnergyAnalyzer"; }

  private:
      // for bonding, we need to create a super-superstructure (cluster)
    xstructure xstr_cluster;                                // cluster of atoms within radius
    std::vector<std::vector<uint>> v_bonded_atom_indices;            // references xstr_cluster
    std::vector<std::vector<uint>> v_nonbonded_atom_indices;         // references xstr_cluster

    bool has_vacancies = false;                                     // are there vacancies present in m_types_config?
    bool bonding_set = false;                                       // have we already found bonding for this configuration?
    double m_energy_uff = AUROSTD_MAX_DOUBLE;

    uint NNDistancesMapPC(uint atom);
    uint NNDistancesMapSC(uint atom);
    void calculateNNDistances(xstructure& xstr, std::vector<uint>& v_vacancies);

    // SERIALIZATION MEMBERS, include macro from POccCalculatorTemplate
#define JSON_POccUFFEnergyAnalyzer_MEMBERS \
  m_initialized, hnf_mat, m_types_config, types2uffparams_map, m_vacancies, m_exploration_radius, distance_matrix, v_dist_nn, xstr_ss, sc2pc_map, pc2sc_map, JSON_POccCalculatorTemplate_MEMBERS
  };
} // namespace pocc

namespace pocc {
  struct GroupMember {
    int basis;  // basis is int in xatom
    aurostd::xvector<double> cpos;
  };
} // namespace pocc

namespace pocc {
  /// @brief Class to ease the construction of @c POccCalculator instances.
  /// Simple builder pattern to currying.
  /// Use the @c with method to chain build arguments and then call @c build to get the constructed @c POccCalculator instance.
  /// For example:
  /// @code POccCalculator pocccalc = POccCalculatorBuilder().with(aflags).with(kflags).with(pocc_flags).build(); @endcode
  ///
  /// @authors
  /// @xlink{https://en.wikipedia.org/wiki/Builder_pattern}
  /// @mod{ST,20250717,created}
  class POccCalculatorBuilder {
  public:
    explicit POccCalculatorBuilder(std::ostream& oss = std::cout) : p_oss(&oss) {}
    explicit POccCalculatorBuilder(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : p_oss(&oss), p_ofs(&FileMESSAGE) {}

    POccCalculatorBuilder& with(const xstructure& xstr_pocc) {
      p_xstr_pocc = &xstr_pocc;
      return *this;
    }
    POccCalculatorBuilder& with(const _aflags& aflags) {
      p_aflags = &aflags;
      return *this;
    }
    POccCalculatorBuilder& with(const _kflags& kflags) {
      p_kflags = &kflags;
      return *this;
    }
    POccCalculatorBuilder& with(const _vflags& vflags) {
      p_vflags = &vflags;
      return *this;
    }
    POccCalculatorBuilder& with(const aurostd::xoption& pocc_flags) {
      p_pocc_flags = &pocc_flags;
      return *this;
    }
    POccCalculatorBuilder& with(std::ostream& oss) {
      p_oss = &oss;
      return *this;
    }
    POccCalculatorBuilder& with(std::ofstream& ofs) {
      p_ofs = &ofs;
      return *this;
    }
    template <typename... BuildArgs> POccCalculatorBuilder& with_all(BuildArgs&&... args) {
      (with(std::forward<BuildArgs>(args)), ...);
      return *this;
    }
    [[nodiscard]] POccCalculator build() const;

    // ST20250715 note: the `with` overloads only work if the arguments are all different types
    // so if you add another xoption to the builder, it either needs a different name or you can inspect the xoption for its specialization somehow

  private:
    std::ostream* p_oss = nullptr;
    std::ofstream* p_ofs = nullptr;
    const _aflags* p_aflags = nullptr;
    const _kflags* p_kflags = nullptr;
    const _vflags* p_vflags = nullptr;
    const xstructure* p_xstr_pocc = nullptr;
    const aurostd::xoption* p_pocc_flags = nullptr;
  };

  class POccCalculator : public POccCalculatorTemplate, public xStream, public JsonSerializable<POccCalculator> {
    friend POccCalculatorBuilder;

  public:
    explicit POccCalculator(std::ostream& oss = std::cout) : xStream(oss) {}
    explicit POccCalculator(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout) : xStream(FileMESSAGE, oss) {}

    template <typename... BuildArgs, aurostd::enable_if_none_of<POccCalculator, BuildArgs...> = true>
    explicit POccCalculator(BuildArgs&&... args) : POccCalculator(build_from(std::forward<BuildArgs>(args)...)) {}
    template <typename... BuildArgs> static POccCalculator build_from(BuildArgs&&... args) { return POccCalculatorBuilder().with_all(std::forward<BuildArgs>(args)...).build(); }

    void clear();

    // inputs
    bool m_initialized = false;
    _kflags m_kflags;                         // standard aflow flags
    _vflags m_vflags;                         // standard aflow flags
    xstructure xstr_sym;                    // will contain symmetry objects (_sym_op, pgroups most important here)
    int n_hnf = 0;
    std::vector<POccUnit> m_pocc_sites;                   // groupings of atoms that are on the same site, non-vacant sites only, relative to xstr_nopocc
    int pocc_atoms_total = 0;
    std::vector<StructureConfiguration> v_str_configs;

    POccUFFEnergyAnalyzer energy_analyzer;
    unsigned long long int hnf_count = 0;
    unsigned long long int types_config_permutations_count = 0;
    unsigned long long int total_permutations_count = 0;
    std::list<POccSuperCellSet> l_supercell_sets;
    // standard flags - ALL options will be handled via xoptions
    aurostd::xoption enumerator_mode;       // how do we determine duplicates - UFF, SNF, ...

    // post-processing
    bool m_convolution = false;
    bool m_count_unique_fast = false;
    std::vector<std::string> m_ARUN_directories;
    double m_Hmix = AUROSTD_MAX_DOUBLE;
    double m_efa = AUROSTD_MAX_DOUBLE;
    int m_temperature_precision = TEMPERATURE_PRECISION;
    int m_zero_padding_temperature = 0;
    bool m_temperatures_int = false;
    uint m_relaxation_max = AUROSTD_MAX_UINT;
    double m_energy_dft_ground = AUROSTD_MAX_DOUBLE;
    uint m_ARUN_directory_ground = AUROSTD_MAX_UINT;
    aurostd::xmatrix<double> m_rdf_all;
    double m_rdf_rmax = AUROSTD_MAX_DOUBLE;
    int m_rdf_nbins = AUROSTD_MAX_INT;
    xDOSCAR m_xdoscar;
    std::vector<double> m_Egap_DOS, m_Egap;
    double m_Egap_DOS_net = AUROSTD_MAX_DOUBLE;
    double m_Egap_net = AUROSTD_MAX_DOUBLE;
    std::vector<std::string> m_vfilenames_xout;  // dielectric file names
    xOUTCAR m_vxout_avg; //averaged dielectric data

    // external methods
    void setPOccFlags(const aurostd::xoption& pocc_flags) override;                    // input flags, e.g., vpflow
    void loadFromAFlags();                                                    // grabs from m_aflags
    void loadFromAFlags(const aurostd::xoption& loader);                      // grabs from m_aflags
    void setPOccStructure(const xstructure& xstr_pocc) override;
    void setAFlags(const _aflags& Aflags) override;                                    // standard _aflags
    void setKFlags(const _kflags& Kflags);                                    // standard _kflags
    void setVFlags(const _vflags& Vflags);                                    // standard _vflags

    // path getters
    std::string getARUNDirectoryPath(uint isupercell) const;
    std::string getOutputPath() const;

    // Modules
    bool inputFilesFoundAnywhereAPL();  // ME20211004
    void createModuleAflowIns(const _xvasp& xvasp, const std::string& MODULE);  // ME20211004

    void writePARTCAR() const;
    void generateStructures(const _xvasp& xvasp);
    xstructure createXStructure(const POccSuperCell& psc, int n_hnf = 0, unsigned long long int hnf_count = 0, unsigned long long int types_config_permutations_count = 0, bool clean_structure = false, bool primitivize = false);
    bool areEquivalentStructuresByUFF(std::list<POccSuperCellSet>::iterator it, const POccSuperCell& psc) const;
    void add2DerivativeStructuresList(const POccSuperCell& psc, std::list<POccSuperCellSet>::iterator i_start, std::list<POccSuperCellSet>::iterator i_end);
    void add2DerivativeStructuresList(const POccSuperCell& psc);
    void getHNFMatSiteConfig(const POccSuperCell& psc, aurostd::xmatrix<double>& hnf_mat, std::vector<std::vector<int>>& v_types_config);
    bool areEquivalentStructures(const POccSuperCell& psc_a, const POccSuperCell& psc_b);
    bool areEquivalentStructures(const xstructure& a, const xstructure& b);
    unsigned long long int runRobustStructureComparison(std::list<POccSuperCellSet>::iterator it);
    void calculateHNF();
    void getTotalPermutationsCount();
    void calculatePOccSuperCellUFF(int thread_id,
                                   std::vector<POccSuperCell>& vpsc,
                                   const std::vector<POccUFFEnergyAnalyzer>& v_energy_analyzer,
                                   const std::vector<std::vector<std::vector<int>>>& vv_types_config,
                                   size_t& npsc_queue,
                                   std::mutex& m_save,
                                   std::mutex& m_job); // SD20230604
    void countUniquePOccSuperCellUFF(int thread_id,
                                     std::map<unsigned long long int, std::unordered_map<long long int, unsigned long int>>& map_unique,
                                     const std::vector<POccSuperCell>& vpsc,
                                     const std::vector<POccUFFEnergyAnalyzer>& v_energy_analyzer,
                                     const std::vector<std::vector<std::vector<int>>>& vv_types_config,
                                     size_t& npsc_queue,
                                     std::mutex& m_save,
                                     std::mutex& m_job); // SD20230609
    void calculate();
    std::string getARUNString(unsigned long long int i) const;
    xstructure getUniqueSuperCell(unsigned long long int i);
    std::vector<xstructure> getUniqueDerivativeStructures();
    std::vector<uint> getMapToPARTCAR(unsigned long long int i, const xstructure& xstr);  // ME20211006
    unsigned long long int getUniqueSuperCellsCount() const;
    // bool printUniqueDerviativeStructures();
    // void resetMaxSLRadius();
    void resetHNFMatrices();
    void resetSiteConfigurations();

    void CleanPostProcessing();
    void convolution();
    void loadDataIntoCalculator();
    void setTemperatureStringParameters();
    void setTemperatureStringParameters(std::vector<double>& v_temperatures);
    void postProcessing();
    void StructuresAllFile2SupercellSets();
    void StructuresUniqueFile2SupercellSets();
    bool QMVASPsFound() const;
    void setDFTEnergies();
    void setEFA();
    void calculateRELAXProperties(double temperature = 300);
    void calculateSTATICProperties(double temperature = 300);
    void calculateDielectricProperties(double temperature = 300);
    void setPOccStructureProbabilities(double temperature = 300); // room temperature
    std::string getTemperatureString(double temperature) const;
    void setAvgRDF(double temperature = 300);  // depends on probabilities
    void setAvgDOSCAR(double temperature = 300);  // depends on probabilities
    void setAvgDielectricData(double temperature = 300);  // depends on probabilities
    void plotAvgDOSCAR(double temperature) const; // no default temperature, needs to be set inside setAvgDOSCAR()
    void plotAvgDOSCAR(const std::string& doscar_path, const std::string& directory = ".") const;
    void plotAvgDOSCAR(const xDOSCAR& xdos, double temperature, const std::string& directory = ".") const;
    void setTempIndependentProperties(aurostd::JSON::object& pocc_out_json) const;
    void setTempDependentProperties(double temperature, aurostd::JSON::object& pocc_out_json) const;

  private:
    static constexpr int A_START = 1;
    static constexpr int C_START = 1;
    static constexpr int F_START = 1;
    static constexpr int B_START = 0;
    static constexpr int D_START = 0;
    static constexpr int E_START = 0;

    // hnf matrices
    int a_start = A_START;
    int b_start = B_START;
    int c_start = C_START;
    int d_start = D_START;
    int e_start = E_START;
    int f_start = F_START;
    std::vector<aurostd::xmatrix<double>> v_unique_superlattices;
    // double max_superlattice_radius;
    // configurations
    aurostd::xmatrix<double> hnf_mat;
    std::vector<std::vector<int>> v_types_config;
    // vector<int> v_config_iterators;
    uint config_iterator = 0;
    std::vector<uint> v_config_order;
    double m_energy_uff_tolerance = DEFAULT_UFF_ENERGY_TOLERANCE;

    void initializePOccStructure();
    void cleanPOccStructure();
    void preparePOccStructure();  // do not write out sym stuff by default

    const aurostd::xmatrix<double>& getLattice() const;
    const std::vector<_sym_op>& getPGroup() const;
    // const StructureConfiguration& getXStrCountConfiguration(uint i) const;

    // void calculateHNF();                   //get n_hnf

    // useful internal methods
    // vector<_sym_op> getPGroups();           //fetch pgroups of xstr_nopocc
    // uint getHNFCount();                     //get count of hnf matrices, refers to v_hnf
    // xmatrix<double> getHNFMatrix(uint i);   //fetch specific hnf matrix, refers to v_hnf
    //_atom getAtom(uint i);                  //grab specific atom, refers to xstr_pocc

    // table stuff
    // set some nice printing precisions and paddings, mostly definitions
    uint getHNFTabelPOCCPrecision() const;
    uint getHNFTableGeneralPrecision() const;
    uint getHNFTableIterationPadding() const;
    uint getHNFTableErrorPadding() const;
    uint getHNFTableColumnPadding() const;
    std::string getHeaderMaxSiteError() const;
    std::string getHeaderMaxStoichError() const;
    std::string getHeaderStoichiometry() const;

    std::string hnfTableHeader();
    void writeHNFTableOutput(int i_hnf, double& stoich_error, double& site_error);
    std::string hnfTableLineOutput(int i_hnf, int str_config);
    // void setHNFTablePadding(int _AFLOW_POCC_PRECISION_);

    void partitionPOccSites();              // get pocc_sites
    aurostd::xvector<double> calculateStoichEachType(std::vector<std::vector<int>>& v_types_config);
    void calculateSymNonPOccStructure(bool verbose = true);          // calculate symmetry of non-pocc structure
    void propagateEquivalentAtoms2POccStructure();  // propagates equivalent atoms info from xstr_sym to xstr_pocc for sorting
    void redecorateXStructures();            // redecorate xstr_nopocc for standardization of UFF energies

    // void getSiteCountConfigurations(int i_hnf,double& stoich_error);
    void getSiteCountConfigurations(int i_hnf);
    void getOptimizedSiteCountConfigurations(int site, int i_hnf, std::vector<POccSiteConfiguration>& v_site_configs);

    // determines site occupancy and vacancy count, given n_hnf
    bool iterateHNFMatrix();                // calculate all unique hnf's
    void setConfigOrder();
    bool getNextSiteConfiguration();
    bool getNextSiteConfiguration(std::vector<std::vector<int>>& v_site_config);

    // CT20200319 - POCC+AEL functions
    void calculateElasticProperties(const std::vector<double>& v_temperatures);
    void setAELOptions(bool& ael_run_postprocess, bool& ael_write_full_results);
    void generateElasticProperties(std::vector<double>& Bvoigt,
                                   std::vector<double>& Breuss,
                                   std::vector<double>& Bvrh,
                                   std::vector<double>& Gvoigt,
                                   std::vector<double>& Greuss,
                                   std::vector<double>& Gvrh,
                                   std::vector<double>& Poisson_ratio,
                                   std::vector<std::vector<std::vector<double>>>& elastic_tensor_list,
                                   std::vector<std::vector<std::vector<double>>>& compliance_tensor_list);
    void getElasticProperties(std::vector<double>& Bvoigt,
                              std::vector<double>& Breuss,
                              std::vector<double>& Bvrh,
                              std::vector<double>& Gvoigt,
                              std::vector<double>& Greuss,
                              std::vector<double>& Gvrh,
                              std::vector<double>& Poisson_ratio,
                              std::vector<std::vector<std::vector<double>>>& elastic_tensor_list,
                              std::vector<std::vector<std::vector<double>>>& compliance_tensor_list);
    void getAverageElasticProperties(const std::vector<double>& v_temperatures,
                                     bool ael_write_full_results,
                                     std::vector<double>& Bvoigt,
                                     std::vector<double>& Breuss,
                                     std::vector<double>& Bvrh,
                                     std::vector<double>& Gvoigt,
                                     std::vector<double>& Greuss,
                                     std::vector<double>& Gvrh,
                                     std::vector<double>& Poisson_ratio,
                                     std::vector<std::vector<std::vector<double>>>& elastic_tensor_list,
                                     std::vector<std::vector<std::vector<double>>>& compliance_tensor_list);

    // CT20200323 - POCC+AGL functions
    void calculateDebyeThermalProperties(const std::vector<double>& v_temperatures);
    void setAGLOptions(bool& agl_run_postprocess, bool& agl_write_full_results); // CT20200722
    void generateDebyeThermalProperties(std::vector<double>& Debye_temperature,
                                        std::vector<double>& Debye_acoustic,
                                        std::vector<double>& Gruneisen,
                                        std::vector<double>& Cv300K,
                                        std::vector<double>& Cp300K,
                                        std::vector<double>& Fvib300K_atom,
                                        std::vector<double>& Fvib300K_cell,
                                        std::vector<double>& Svib300K_atom,
                                        std::vector<double>& Svib300K_cell,
                                        std::vector<double>& kappa300K,
                                        std::vector<std::vector<double>>& agl_temperatures,
                                        std::vector<std::vector<double>>& agl_gibbs_energies_atom,
                                        std::vector<std::vector<double>>& agl_vibrational_energies_atom); // CT20200722
    void getDebyeThermalProperties(std::vector<double>& Debye_temperature,
                                   std::vector<double>& Debye_acoustic,
                                   std::vector<double>& Gruneisen,
                                   std::vector<double>& Cv300K,
                                   std::vector<double>& Cp300K,
                                   std::vector<double>& Fvib300K_atom,
                                   std::vector<double>& Fvib300K_cell,
                                   std::vector<double>& Svib300K_atom,
                                   std::vector<double>& Svib300K_cell,
                                   std::vector<double>& kappa300K,
                                   std::vector<std::vector<double>>& agl_temperatures,
                                   std::vector<std::vector<double>>& agl_gibbs_energies_atom,
                                   std::vector<std::vector<double>>& agl_vibrational_energies_atom);
    void getAverageDebyeThermalProperties(const std::vector<double>& v_temperatures,
                                          bool agl_write_full_results,
                                          std::vector<double>& Debye_temperature,
                                          std::vector<double>& Debye_acoustic,
                                          std::vector<double>& Gruneisen,
                                          std::vector<double>& Cv300K,
                                          std::vector<double>& Cp300K,
                                          std::vector<double>& Fvib300K_atom,
                                          std::vector<double>& Fvib300K_cell,
                                          std::vector<double>& Svib300K_atom,
                                          std::vector<double>& Svib300K_cell,
                                          std::vector<double>& kappa300K,
                                          std::vector<std::vector<double>>& agl_temperatures,
                                          std::vector<std::vector<double>>& agl_gibbs_energies_atom,
                                          std::vector<std::vector<double>>& agl_vibrational_energies_atom);

    // ME20210927 - APL functions
    void calculatePhononPropertiesAPL(const std::vector<double>& v_temperatures);
    std::vector<apl::PhononCalculator> initializePhononCalculators();
    std::vector<xDOSCAR> getPhononDoscars(std::vector<apl::PhononCalculator>& vphcalc, aurostd::xoption& dosopts, std::vector<int>& vexclude);
#ifdef AFLOW_MULTITHREADS_ENABLE
    void calculatePhononDOSThread(uint i, const std::vector<uint>& vcalc, const aurostd::xoption& aplopts, std::vector<apl::DOSCalculator>& vphdos, std::vector<xDOSCAR>& vxdos, std::mutex& m);
#else
    void calculatePhononDOSThread(uint i, const std::vector<uint>& vcalc, const aurostd::xoption& aplopts, std::vector<apl::DOSCalculator>& vphdos, std::vector<xDOSCAR>& vxdos);
#endif
    xDOSCAR getAveragePhononDos(double T, const std::vector<xDOSCAR>& vxdos);
    // AS20210204 QHA
    void calculateQHAProperties();
    void calculateQHAPropertiesAVG(const std::vector<double>& v_temperatures);

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccCalculator deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccCalculator"; }

    // SERIALIZATION MEMBERS Avoiding: xoption enumerator_mode, include macro from POccCalculatorTemplate
#define JSON_POccCalculator_MEMBERS                                                                                                                                                                           \
  m_initialized, m_kflags, m_vflags, xstr_sym, n_hnf, m_pocc_sites, pocc_atoms_total, energy_analyzer, hnf_count, types_config_permutations_count, v_str_configs, total_permutations_count, l_supercell_sets, \
      m_convolution, m_count_unique_fast, m_ARUN_directories, m_Hmix, m_efa, m_temperature_precision, m_zero_padding_temperature, m_temperatures_int, m_relaxation_max, m_energy_dft_ground,                  \
      m_ARUN_directory_ground, m_rdf_all, m_rdf_rmax, m_rdf_nbins, m_xdoscar, m_Egap_DOS, m_Egap_DOS_net, m_vfilenames_xout, a_start, b_start, c_start, d_start, e_start, f_start,  \
      v_unique_superlattices, hnf_mat, v_types_config, config_iterator, v_config_order, m_energy_uff_tolerance, JSON_POccCalculatorTemplate_MEMBERS
  };
} // namespace pocc

namespace pocc {
  class POccStructuresFile : public xStream, public JsonSerializable<POccStructuresFile> {
  public:
    POccStructuresFile(std::ostream& oss = std::cout);
    POccStructuresFile(std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    POccStructuresFile(const std::string& fileIN, std::ostream& oss = std::cout);
    POccStructuresFile(const std::string& fileIN, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    POccStructuresFile(const _aflags& aflags, std::ostream& oss = std::cout);
    POccStructuresFile(const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);
    POccStructuresFile(const std::string& fileIN, const _aflags& aflags, std::ostream& oss = std::cout);
    POccStructuresFile(const std::string& fileIN, const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss = std::cout);

    void clear();

    bool m_initialized = false;
    std::string m_content;
    std::vector<std::string> m_vcontent;
    std::string m_filename;
    _aflags m_aflags;
      // vector<string> m_ARUN_directories;
    std::list<POccSuperCellSet> l_supercell_sets;
    aurostd::xoption m_fileoptions;
    std::vector<std::vector<std::vector<uint>>> m_vPOSCAR_lines;  // contains start/stop indices for POSCARs in m_vcontent

      // initializers
    bool initialize(std::ostream& oss);
    bool initialize(std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize();
    bool initialize(const std::string& fileIN, std::ostream& oss);
    bool initialize(const std::string& fileIN, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const std::string& fileIN);
    bool initialize(const _aflags& aflags, std::ostream& oss);
    bool initialize(const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const _aflags& aflags);
    bool initialize(const std::string& fileIN, const _aflags& aflags, std::ostream& oss);
    bool initialize(const std::string& fileIN, const _aflags& aflags, std::ofstream& FileMESSAGE, std::ostream& oss);
    bool initialize(const std::string& fileIN, const _aflags& aflags);

    void setAFlags(const _aflags& aflags);
    void readFile(const std::string& fileIN);
    void processFile();
    bool getARUNDirectories(std::vector<std::string>& ARUN_directories, bool tryDirectoryLS = true);
    bool loadDataIntoCalculator(POccCalculator& pcalc, bool tryDirectoryLS = true);

    friend std::ostream& operator<<(std::ostream&, const POccStructuresFile&);

    // JSON load/dump
  protected:
    [[nodiscard]] aurostd::JSON::object serialize() const override;
    POccStructuresFile deserialize(const aurostd::JSON::object& jo) override;
    [[nodiscard]] std::string getJsonID() const override { return "POccStructuresFile"; }

  private:
    void free();

    // SERIALIZATION MEMBERS, not serializaing xoption vars
#define JSON_POccStructuresFile_MEMBERS m_initialized, m_content, m_aflags, l_supercell_sets, m_vPOSCAR_lines
  };
} // namespace pocc

namespace pocc {
  /// This class is used to calculate  thermodynamic properties over the ensemble,
  /// defined by structures generated using POCC method.
  class EnsembleThermo : public xStream {
  public:
    EnsembleThermo(std::ostream& oss = std::cout);
    EnsembleThermo(const std::string& currentDir,
                   std::vector<std::string>& directories,
                   const std::string& filename,
                   const std::string& calc_type,
                   apl::EOSmethod eos_method,
                   bool isFVTprovided,
                   std::ofstream& FileMESSAGE,
                   std::ostream& oss = std::cout);
    apl::QHA qha;
    apl::EOSmethod eos_method = apl::EOS_SJ;
    uint Nstructures = 0;
    int Nvolumes = 0;
    int nrows = 0;
    double Ensemble_Vmin = AUROSTD_MAX_DOUBLE;
    double Ensemble_Vmax = 0;
    std::string currentDirectory;
    aurostd::xvector<double> T;
    aurostd::xmatrix<double> FV;
    aurostd::xvector<double> volumes;
    std::vector<int> degeneracies;
    std::vector<aurostd::xmatrix<double>> coeffs_list;
    aurostd::xvector<double> Veq, Feq, B, Bprime, Cv, Cp, gamma, beta;
    double logZ(const aurostd::xvector<double>& E, const std::vector<int>& degeneracies, double T);
    aurostd::xvector<double> calcThermalExpansionSG(const aurostd::xvector<double>& volumes, double dT);
    aurostd::xvector<double> calcIsobaricSpecificHeatSG(const aurostd::xvector<double>& free_energies, double dT);
    void calculateThermodynamicProperties();
    void writeThermodynamicProperties() const;
    void clear();

  private:
    void readFVTParameters(const std::string& filename, const std::string& blockname, uint& Nvolumes, uint& Ntemperatures);
    void readFVTdata(const std::string& dirname, const std::string& filename, const std::string& blockname, uint n_volumes, uint n_temperatures, aurostd::xvector<double>& t, aurostd::xmatrix<double>& c, double& Vmin, double& Vmax);
    bool readCoeffData(const std::string& filename, const std::string& blockname, aurostd::xvector<double>& T, aurostd::xmatrix<double>& coeffs);
    void readCoeffParameters(const std::string& filename, double& Vmin, double& Vmax);

    void free();
  };
} // namespace pocc

#endif  // _AFLOW_POCC_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
