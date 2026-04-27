
#ifndef AFLOW_XSERIALIZATION_TPP
#define AFLOW_XSERIALIZATION_TPP

#include "config.h"

#include <string>

#include "aurostd.h"
#include "aurostd_time.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"
#include "aurostd_xparser_json.h"

/// @return The class serialized as a Json object
/// @authors
/// @mod{ST,20250108,created}
template <typename T> [[nodiscard]] aurostd::JSON::object JsonSerializable<T>::dumpToJson() const {
  aurostd::JSON::object result = this->serialize();
  result.join(dumpJsonMeta());
  return result;
}

/// @return The class serialized as a string representation of a Json object
/// @authors
/// @mod{ST,20250108,created}
template <typename T> [[nodiscard]] std::string JsonSerializable<T>::dumpToString() const {
  return dumpToJson().toString();
}

/// @brief Writes the class to file, serialized as a string representation of a Json object
/// @param path the path to write to
/// @param compression the compression mode to use, may use aurostd::compression_type::None or 0 for no compression
/// @authors
/// @mod{ST,20250108,created}
template <typename T> void JsonSerializable<T>::dumpToFile(const std::filesystem::path& path, aurostd::compression_type compression) const {
  dumpToJson().saveFile(path, compression);
}

/// @param jo The json object to load from
/// @return deserailized instance
/// @authors
/// @mod{ST,20250108,created}
template <typename T> T JsonSerializable<T>::loadFromJson(const aurostd::JSON::object& jo) {
  T t;
  JsonSerializable& tr = t;
  tr.checkClassID(jo);
  return tr.deserialize(jo);
}

/// @param path The path to read the string representation of json object to load from
/// @return deserailized instance
/// @authors
/// @mod{ST,20250108,created}
template <typename T> T JsonSerializable<T>::loadFromFile(const std::filesystem::path& path) {
  return loadFromJson(aurostd::JSON::loadFile(path));
}

/// @param content The string representation of json object to load from
/// @return deserailized instance
/// @authors
/// @mod{ST,20250108,created}
template <typename T> T JsonSerializable<T>::loadFromString(const std::string& content) {
  return loadFromJson(aurostd::JSON::loadString(content));
}

/// @return The metadata json object for json serialization
/// @authors
/// @mod{ST,20250108,created}
template <typename T> [[nodiscard]] aurostd::JSON::object JsonSerializable<T>::dumpJsonMeta() const {
  return aurostd::JSON::object(aurostd::JSON::Dictionary{
      {"_aflow_serialization", aurostd::JSON::object(aurostd::JSON::Dictionary{
                                   {"datetime", aurostd::get_datetime_formatted("-", false)},
                                   {"aflow_version", AFLOW_VERSION},
                                   {"class_id", getJsonID()},
                               })}
  });
}

/// @brief Checks the json object to see if its class id matches that of the specialized caller.
/// Throws if no match, does nothing if matches.
/// @throw aurostd::xerror when the class id does not match
/// @authors
/// @mod{ST,20250108,created}
template <typename T> void JsonSerializable<T>::checkClassID(const aurostd::JSON::object& jo) const {
  if (const std::string id(jo["_aflow_serialization"]["class_id"]); id != getJsonID()) {
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Cannot load mismatched class id. Expected: " + getJsonID() + ", Got: " + id);
  }
}

#endif // AFLOW_XSERIALIZATION_TPP
