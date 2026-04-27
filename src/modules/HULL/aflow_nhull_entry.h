#ifndef AFLOW_ENTRY_H
#define AFLOW_ENTRY_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xhttp.h"
#include "AUROSTD/aurostd_xserialization.h"
#include "AUROSTD/aurostd_xvector.h"

#include "aflowlib/aflowlib_web_interface.h"

namespace nhull {
  //struct represents all error flags that can exclude entries from calculation
  //All global entry flags are updated in filterEntries() in nhull.cpp
  struct Flags : public JsonSerializable<Flags> {
    bool extraneous_energy_entry; // indicates whether entry was flagged for having anomalous energy value
    bool incompatible_entry; // indicates whether entry was loaded with species not corresponding to input vector (this would be an error)
    bool plus_U_used;
    bool outside_IQR; // indicates that entry was flagged for having enthalpy outside IQR range filter in this run
    bool duplicate_entry;
    bool vspecie_vcomposition_size_error; //true if vspecie.size() != vcomp.size()
    bool pocc_prototype; //true if this entry is a pocc prototype
    bool uses_disallowed_dft_type;
    bool uses_NUPDOWN;
    bool arun_entry;
    bool bad_database;
    bool corrected_energy; //indicates that enthalpies were manually changed. DO NOT THROW FLAG FOR THIS TYPE IN is_flagged()
    bool excluded_stability_criterion;
    bool excluded_nminus1;

    //default constructor- initialize all flags to false
    Flags() :
        extraneous_energy_entry(false),
        incompatible_entry(false),
        plus_U_used(false),
        outside_IQR(false),
        duplicate_entry(false),
        vspecie_vcomposition_size_error(false),
        pocc_prototype(false),
        uses_disallowed_dft_type(false),
        uses_NUPDOWN(false),
        arun_entry(false),
        bad_database(false),
        corrected_energy(false),
        excluded_nminus1(false),
        excluded_stability_criterion(false) {}

    //returns true if any flag is true
    //if pocc_run = true, ignores filter on pocc compounds
    [[nodiscard]] bool is_flagged(const bool pocc_run) const {
      bool flagged = false;
      flagged |= this->extraneous_energy_entry;
      flagged |= this->incompatible_entry;
      flagged |= this->plus_U_used;
      flagged |= this->outside_IQR;
      flagged |= this->duplicate_entry;
      flagged |= this->vspecie_vcomposition_size_error;
      if (!pocc_run) {
        flagged |= this->pocc_prototype;
      }
      flagged |= this->uses_disallowed_dft_type;
      flagged |= this->uses_NUPDOWN;
      flagged |= this->arun_entry;
      flagged |= this->bad_database;
      flagged |= this->excluded_stability_criterion;
      flagged |= this->excluded_nminus1;
      return flagged;
    }

  protected:
    [[nodiscard]] std::string getJsonID() const override { return "Flags"; };

  public:
    ~Flags() override = default;

    // JSON load/dump
    inline static const std::string JSONID = "Flags";

    [[nodiscard]] aurostd::JSON::object serialize() const override;

    Flags deserialize(const aurostd::JSON::object& jo) override;

    // SERIALIZATION MEMBERS
#define JSON_Flags_MEMBERS extraneous_energy_entry, incompatible_entry, plus_U_used, outside_IQR, duplicate_entry
  };

  // general object for a single alloy
  class Entry : public JsonSerializable<Entry> {
  public:
    bool m_initialized;
    aurostd::xvector<double> m_coord;           // most general hull coordinates full coordinate (non-stoichiometric coord with first index removed)
    bool m_is_artificial;
    bool m_has_entry;  //indexing variable
    uint m_i_coord_group;
    uint nh_index;   //indexing variable for performance and ease of calling, idexed in order of loading from aflowlib
    double m_stability_criterion; //distance to hull with a chosen auid removed from hull
    double m_nminus1_enthalpy_gain; //distance to hull with all species below dimension n removed
    double relative_stability_criterion; //stability criterion divided by enthalpy_formation_atom
    std::string aflow_lib_entry;
    aflowlib::_aflowlib_entry lib_entry_direct;
    std::string compound;
    std::vector<double> vcomposition;
    std::vector<std::string> vspecies;//specie vector for just this entry
    std::vector<std::string> global_vspecies; //specie vector of the full dimension hull run
    std::string prototype;
    std::string auid;
    std::string aurl;
    int nspecies;
    int pocc_degeneracy; // DG value from a pocc entry

    double spin_atom;
    double enthalpy_formation_atom;
    double entropic_temperature;

    double distance_hull_enthalpy_formation_atom;
    Flags flags; //all flags associated with entry
    bool is_hull_point; // indicates whether point is directly on hull
    bool lowest_energy_alloy; // indicates whether point is lowest energy of its reduced stoich family
    bool flagged_entry; // indicates that entry failed one of the filtering steps if true
    bool m_uses_cce; //true if entry enthalpy corrected with CCE
    std::map<std::string, double> nhull_phase_decomp; //used as a convenient storage for JSON export
    std::map<Entry, double> m_nhull_phase_decomp_entries;
    std::map<std::string, int> compound_map; // map representation of compound string, Ex: "Mn2Pd4" -> {{"Mn", 1}, {"Pd", 2}}; this is in reduced stoichiometry
    std::string ldau_TLUJ; // plus U calculation value

    static constexpr double EXTRANEOUS_ENERGY_LIMIT = 1000;

    // constructors
    Entry() = default;
    Entry(const aurostd::xvector<double>& coord, bool is_artificial);
    Entry(const aflowlib::_aflowlib_entry& aflow_lib_entry, const std::vector<std::string>& global_specie_vector, int nhull_index, bool neglect_cce);
    Entry(const aurostd::xvector<double>& coord, std::ostream& oss = std::cout, bool has_stoich_coords = false, bool formation_energy_coord = false, bool is_artificial = false, const uint ch_index = 0); //constructor used only for GFA:
    Entry(const std::map<std::string, double>& nhull_phase_decomp, const std::string& auid, const double dhull);

    // operators

    bool operator==(const Entry& other_entry) const { return auid == other_entry.auid; }
    bool operator!=(const Entry& other_entry) const { return auid != other_entry.auid; }

    bool operator<(const Entry& other_entry) const { return auid < other_entry.auid; } //necessary for map initialization (weak ordering)

    // constructor utilities

    void compoundMapInit();

    // methods about points

    static std::vector<Entry> pointsToEntries(const std::vector<aurostd::xvector<double>>& points);
    [[nodiscard]] static bool getEntryFromPoint(const std::vector<Entry>& entries, Entry& found_entry, const aurostd::xvector<double>& point);
    [[nodiscard]] bool entryEqualsPoint(const aurostd::xvector<double>& point) const;
    static bool isPointUnary(const aurostd::xvector<double>& point);
    [[nodiscard]] aurostd::xvector<double> compoundToPoint(const std::vector<std::string>& sorted_name_store, bool full_point) const;

    void init_mCoord(std::vector<std::string> specie_vector);
    [[nodiscard]] static std::vector<aurostd::xvector<double>> entriesToPoints(std::vector<Entry>);

    // utilities

    [[nodiscard]] std::vector<std::string> toSpecieVector() const; // creates representation of compound string, Ex: "Mn1Pd2" -> {{"Mn", 1}, {"Pd", 2}}; this is run in constructor
    void enthalpyUnitConversion();
    static uint maxDimension(const std::vector<Entry>& entries);
    static std::map<std::string, int> compoundMap(std::string compound);
    static double minEnergy(const std::vector<Entry>& entries);
    bool plusUCheck();
    bool isCompoundCompatible(const std::vector<std::string>& alloy_name_store);
    bool isCompoundCompatible(const aurostd::xvector<int>& alloy_index_store);
    static bool findAUID(const std::string auid, const std::vector<Entry>& entries);
    [[nodiscard]] static std::vector<Entry> getLowestEnergyEntries(const std::vector<Entry>& entries, int filter_out_below_dim = 0);
    [[nodiscard]] static std::vector<Entry> getLowestEnergyNminus1Entries(const std::vector<Entry>& entries, int filter_out_above_dim);
    static std::vector<Entry> getLowestEnergySCEntries(const std::vector<Entry>& entries, const Entry& filter_out_entry);

    // setters

    void set_distance_hull_enthalpy_formation_atom(double distance_hull_enthalpy_formation_atom);
    void set_relative_stability_criterion(double distance_hull_enthalpy_formation_atom, double enthalpy_formation_atom);

    void set_pocc_degeneracy(int pocc_degeneracy) { this->pocc_degeneracy = pocc_degeneracy; }

    // getters

    [[nodiscard]] bool uses_LDAU() const { return !ldau_TLUJ.empty(); }
    [[nodiscard]] bool is_outside_IRQ(const double inner_quartile_cutoff);
    [[nodiscard]] bool is_energy_extraneous() const { return enthalpy_formation_atom > EXTRANEOUS_ENERGY_LIMIT; }
    [[nodiscard]] std::map<std::string, double> phaseDecompToMap() const;
    [[nodiscard]] static Entry getEntryFromAuid(const std::vector<Entry>& entries, const std::string auid);
    [[nodiscard]] aurostd::xvector<double> getStoichiometricCoords() const;
    [[nodiscard]] double getLastCoord() const;
    [[nodiscard]] uint getCoordDimension() const;
    [[nodiscard]] aurostd::xvector<double> findPhaseDecompRatios() const;
    [[nodiscard]] aurostd::xvector<int> getElementsPresent() const;
    [[nodiscard]] static aurostd::xvector<double> getStoichCoords(aurostd::xvector<double> m_coord);

    [[nodiscard]] bool isUnary() const;
    [[nodiscard]] bool isWithinHalfHull(bool lower_hull) const;
    [[nodiscard]] bool isGState() const;
    [[nodiscard]] bool entryIdentical(const Entry& other) const;
    static bool entryUnique(const std::vector<uint>& unique_entries, Entry& other, const std::vector<Entry>& entries);

  protected:
    [[nodiscard]] std::string getJsonID() const override { return "Entry"; };

  public:
    ~Entry() override = default;

    // JSON load/dump
    inline static const std::string JSONID = "Entry";

    [[nodiscard]] aurostd::JSON::object serialize() const override;

    Entry deserialize(const aurostd::JSON::object& jo) override;

    // SERIALIZATION MEMBERS
#define JSON_Entry_MEMBERS                                                                                                                                                                                     \
  compound, flagged_entry, prototype, auid, nspecies, spin_atom, enthalpy_formation_atom, entropic_temperature, lowest_energy_alloy, distance_hull_enthalpy_formation_atom, is_hull_point, nhull_phase_decomp, \
      flags, m_stability_criterion, m_nminus1_enthalpy_gain, relative_stability_criterion
  };

  std::vector<aurostd::xvector<double>> entriesToPoints(const std::vector<Entry>& entries, const std::vector<std::string>& specie_vector);

  std::vector<std::vector<double>> entriesToPoints(const std::map<std::map<std::string, int>, double>& entry_map, const std::vector<std::string>& specie_vector);
} // namespace nhull

#endif // AFLOW_ENTRY_H
