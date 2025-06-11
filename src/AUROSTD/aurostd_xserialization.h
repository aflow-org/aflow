// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                                                                         *
// ***************************************************************************

#ifndef AFLOW_XSERIALIZATION_H
#define AFLOW_XSERIALIZATION_H

#include <filesystem>
#include <string>

#include "aurostd_xfile.h"

namespace aurostd::JSON {
  struct object;
}

/// @brief Interface for Json serialization functionality using CRTP.
/// @tparam T The serializable type and inheriting class
///
/// Inheriting this class provides json serialization capabilities for reading/writing json objects/strings/files.
/// This class should be inherited with public inheritance in order for it to be of use.
/// This class uses CRTP (Curiously Recurring Template Pattern) @link{https://en.cppreference.com/w/cpp/language/crtp}
/// in order to have static methods which act as constructors for the derived (inheriting) class.
/// Inheriting this class will require the override implementation of three protected methods to serialize,
/// deserialize, and check identity.
///
/// Requirements on the templated class T:
/// 1. Must be default constructible.
/// 2. Must inherit this class.
///
/// The implemntation of serialize/deserialize should write/read the members you wish to serialize to/from
/// aurostd::JSON:object. Implementation of the third ID method should return a string ID for the class, typically
/// the class name. See @ref MockSerializableMember and @ref MockSerializableMain for example implementations.
///
/// Note: It may be useful to use a macro to define a comma separated list of the members the class will serialize.
/// If defined, this macro should have the name JSON_<class ID/class name>_MEMBERS. This macro would then be passed to
/// the `AST_JSON_SETTER` and `AST_JSON_GETTER` macros to generate assignments for deserialization and a dictionary/map
/// for serialization, respectively. Again, see @ref MockSerializableMain for an example. These macros will only work
/// if the type is supported for cast conversion by @ref aurostd::JSON::object class, which includes many numeric,
/// list-like, and dictionary-like types as well as others inheriting JsonSerializable.
///
/// @authors
/// @mod{ST,20250108,created}
/// @mod{ST,20250315,slightly modify interface for clearer contract on inheriters}
template <typename T /* <: */> class JsonSerializable {
  friend aurostd::JSON::object;

public:
  virtual ~JsonSerializable() = default;

  [[nodiscard]] aurostd::JSON::object dumpToJson() const;
  [[nodiscard]] std::string dumpToString() const;
  void dumpToFile(const std::filesystem::path &path, aurostd::compression_type compression = aurostd::compression_type::None) const;

  static T loadFromJson(const aurostd::JSON::object &jo);
  static T loadFromFile(const std::filesystem::path &path);
  static T loadFromString(const std::string &content);

protected:
    /// Performs the actual deserialization of members
  virtual T deserialize(const aurostd::JSON::object &jo) = 0;
    /// Performs the actual serialization of members
  [[nodiscard]] virtual aurostd::JSON::object serialize() const = 0;
    /// Returns the string ID of the class
  [[nodiscard]] virtual std::string getJsonID() const = 0;

private:
  [[nodiscard]] aurostd::JSON::object dumpJsonMeta() const;
  void checkClassID(const aurostd::JSON::object &jo) const;
};

#include "aurostd_xserialization.tpp" // NOLINT / template implementation

#endif // AFLOW_XSERIALIZATION_H
