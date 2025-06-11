// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// This JSON class is the evolution of different prior solutions to integrate JSON with the AFLOW source base.
// hagen.eckert@duke.edu

#ifndef _AUROSTD_XPARSER_JSON_CPP_
#define _AUROSTD_XPARSER_JSON_CPP_

#include "aurostd_xparser_json.h"

#include <cmath>
#include <codecvt>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "aurostd.h"
#include "aurostd_xerror.h"
#include "aurostd_xfile.h"
#include "aurostd_xscalar.h"

namespace aurostd {

  /// @struct JSON::object
  /// @brief storge container for a JSON object
  ///
  /// @authors
  /// @mod{HE,20220924,created struct}
  ///
  /// @see
  /// @xlink{"JSON definition",https://www.json.org/json-en.html}
  /// @xlink{"Parsing JSON is a Minefield",https://seriot.ch/projects/parsing_json.html}

  /// @brief direct index access to JSON::object_types::LIST objects
  /// @param index list index
  /// @return JSON::object stored at index
  /// @note throws an error if used on something other than JSON::object_types::LIST
  /// @note the returned JSON::object can be directly be saved into an appropriate variable
  /// like `string content = json_obj[3];`
  ///
  /// @authors
  /// @mod{HE,20220924,created}
  JSON::object &JSON::object::operator[](const size_t index) const {
    if (this->type == JSON::object_types::LIST) {
      const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(this->obj);
      return content->operator[](index);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a list", _INDEX_ILLEGAL_);
    }
  }

  /// @brief direct key access to JSON::object_types::DICTIONARY objects
  /// @param key dictionary key
  /// @return JSON::object stored at key
  /// @note throws an error if used on something other than JSON::object_types::DICTIONARY
  /// @note the returned JSON::object can be directly be saved into an appropriate variable
  /// like `string content = json_obj['my_key'];`
  ///
  /// @authors
  /// @mod{HE,20220924,created}
  JSON::object &JSON::object::operator[](const std::string &key) const {
    if (this->type == JSON::object_types::DICTIONARY) {
      const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(this->obj);
      return content->operator[](key);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a dictionary", _INDEX_ILLEGAL_);
    }
  }

  /// @brief insertion operator: for creating toString
  /// @note does not work externally - the implicit conversion of jo would interfere with the rest of the code
  ostream &operator<<(ostream &os, const JSON::object &jo) {
    os << jo.toString();
    return os;
  }

  /// @brief direct key access to JSON::object_types::DICTIONARY objects
  JSON::object &JSON::object::operator[](const char *key) const {
    if (this->type == JSON::object_types::DICTIONARY) {
      const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(this->obj);
      return content->operator[](key);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "Not a dictionary", _INDEX_ILLEGAL_);
    }
  }

  /// @brief converts a JSON::object into a string
  /// @param json_format if `true` add encapsulate strings in `"` (default: `true`)
  /// @param escape_unicode if `true` escape unicode in strings - `αβγ` becomes `\u03b1\u03b2\u03b3` (default: `true`)
  /// @return JSON formatted string
  ///
  /// @authors
  /// @mod{HE,20220924,created}
  std::string JSON::object::toString(const bool json_format, const bool escape_unicode) const {
    bool first = true;
    stringstream result;
    switch (type) {
      case object_types::DICTIONARY: {
        result << "{";
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        for (const auto &entry : *content) {
          if (first) {
            first = false;
          } else {
            result << ",";
          }
          result << "\"" << entry.first << "\":" << entry.second; //
        }
        result << "}";
      } break;
      case object_types::LIST: {
        result << "[";
        const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
        for (const auto &entry : *content) {
          if (first) {
            first = false;
          } else {
            result << ",";
          }
          result << entry; //
        }
        result << "]";
      } break;
      case object_types::STRING: {
        const std::shared_ptr<std::string> content = std::static_pointer_cast<std::string>(obj);
        if (json_format) {
          result << "\"" << JSON::escape(*content, escape_unicode) << "\"";
        } else {
          return *content;
        }
      } break;
      case object_types::INTEGER: {
        const std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
        result << *content;
      } break;
      case object_types::FLOAT: {
        const std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
        result << *content;
      } break;
      case object_types::T: {
        result << "true";
      } break;
      case object_types::F: {
        result << "false";
      } break;
      case object_types::NONE: {
        result << "null";
      } break;
      default: break;
    }
    return result.str();
  }

  /// @brief save JSON::object to file
  /// @param file_path path to save to
  /// @param escape_unicode if `true` escape unicode in strings - `αβγ` becomes `\u03b1\u03b2\u03b3` (default: `true`)
  /// @param ct compression type
  /// @authors
  /// @mod{HE,20221110,created}
  void JSON::object::saveFile(const std::string &file_path, const compression_type ct, const bool escape_unicode) const {
    aurostd::string2file(this->toString(true, escape_unicode), file_path, ct);
  }

  /// @brief change this JSON::object to a JSON::object_types::STRING
  /// @authors
  /// @mod{HE,20220926,created}
  void JSON::object::fromString(const std::string &content) {
    const std::shared_ptr<std::string> new_string = std::make_shared<std::string>();
    *new_string = content;
    this->obj = new_string;
    this->type = object_types::STRING;
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST
  /// @authors
  /// @mod{ST,20241218,created}
  void JSON::object::fromList(const std::vector<JSON::object> &content) {
    const std::shared_ptr<JSON::List> new_list = std::make_shared<JSON::List>();
    *new_list = content;
    this->obj = new_list;
    this->type = object_types::LIST;
  }

  /// @brief change this JSON::object to a JSON::object_types::DICTIONARY
  /// @authors
  /// @mod{ST,20241218,created}
  void JSON::object::fromDictionary(const std::map<std::string, JSON::object> &content) {
    const std::shared_ptr<JSON::Dictionary> new_dict = std::make_shared<JSON::Dictionary>();
    *new_dict = content;
    this->obj = new_dict;
    this->type = object_types::DICTIONARY;
  }

  /// @brief converting constructor: set the content of this JSON::object based on a char
  JSON::object::object(const char *content) {
    this->fromString(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a string
  JSON::object::object(const std::string &content) {
    this->fromString(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a bool
  JSON::object::object(bool content) {
    if (content) {
      this->type = object_types::T;
    } else {
      this->type = object_types::F;
    }
  }

  /// @brief converting constructor: set the content of this JSON::object to null
  /// @note pass in a nullptr
  JSON::object::object(std::nullptr_t content) {
    this->type = object_types::NONE;
    this->obj = content;
  }

  /// @brief create an empty JSON::object of a given JSON::object_types
  /// @authors
  /// @mod{HE,20221031,created}
  JSON::object::object(object_types create_type) {
    this->type = create_type;
    switch (create_type) {
      {
        case object_types::DICTIONARY: {
          const std::shared_ptr<JSON::Dictionary> content = std::make_shared<JSON::Dictionary>();
          this->obj = content;
          break;
        }
        case object_types::LIST: {
          const std::shared_ptr<JSON::List> content = std::make_shared<JSON::List>();
          this->obj = content;
          break;
        }
        case object_types::STRING: {
          const std::shared_ptr<std::string> content = std::make_shared<std::string>();
          this->obj = content;
          break;
        }
        case object_types::INTEGER: {
          const std::shared_ptr<long long int> content = std::make_shared<long long int>();
          this->obj = content;
          break;
        }
        case object_types::FLOAT: {
          const std::shared_ptr<double> content = std::make_shared<double>();
          this->obj = content;
          break;
        }
        default: {
          this->obj = nullptr;
        }
      }
    }
  }

  JSON::object::object(const std::vector<JSON::object> &content) {
    this->fromList(content);
  }

  JSON::object::object(const std::map<std::string, JSON::object> &content) {
    this->fromMap(content);
  }

  ///@brief assignment operator for char
  JSON::object &JSON::object::operator=(const char *content) {
    this->fromString(content);
    return *this;
  }

  ///@brief assignment operator for string
  JSON::object &JSON::object::operator=(const std::string &content) {
    this->fromString(content);
    return *this;
  }

  ///@brief assignment operator for bool
  JSON::object &JSON::object::operator=(bool content) {
    if (content) {
      this->type = object_types::T;
    } else {
      this->type = object_types::F;
    }
    return *this;
  }

  ///@brief assignment operator for nullptr
  JSON::object &JSON::object::operator=(std::nullptr_t content) {
    this->type = object_types::NONE;
    this->obj = content;
    return *this;
  }

  ///@brief assignment operator for JSON::List
  JSON::object &JSON::object::operator=(const List &content) {
    this->fromList(content);
    return *this;
  }

  ///@brief assignment operator for JSON::Dictionary
  JSON::object &JSON::object::operator=(const Dictionary &content) {
    this->fromDictionary(content);
    return *this;
  }

  ///@brief conversion function for bool
  ///@note for LIST, DICTIONARY and STRING tests emptiness
  ///@note for NULL returns always false
  JSON::object::operator bool() const {
    switch (type) {
      {
        case object_types::DICTIONARY: {
          const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
          if (content->empty()) {
            return false;
          } else {
            return true;
          }
        }
        case object_types::LIST: {
          const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
          if (content->empty()) {
            return false;
          } else {
            return true;
          }
        }
        case object_types::STRING: {
          const std::shared_ptr<std::string> content = std::static_pointer_cast<std::string>(obj);
          if (content->empty()) {
            return false;
          } else {
            return true;
          }
        }
        case object_types::INTEGER: {
          const std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          if (*content) {
            return true;
          } else {
            return false;
          }
        }
        case object_types::FLOAT: {
          const std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          if (static_cast<bool>(*content)) {
            return true;
          } else {
            return false;
          }
        }
        case object_types::T: {
          return true;
        }
        case object_types::F: {
          return false;
        }
        case object_types::NONE: {
          return false;
        }
        default: {
          return false;
        }
      }
    }
  }

  ///@brief conversion function for double
  JSON::object::operator double() const {
    switch (type) {
      {
        case object_types::INTEGER: {
          const std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          return (double) *content;
        }
        case object_types::FLOAT: {
          const std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          return *content;
        }
        case object_types::T: return 1.0;
        case object_types::F: return 0.0;
        case object_types::NONE: return NAN;
        default: throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON double conversion failed: is not a number: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief conversion function for float
  JSON::object::operator long double() const {
    return static_cast<double>(*this);
  }

  ///@brief conversion function for float
  JSON::object::operator float() const {
    return static_cast<double>(*this);
  }

  ///@brief conversion function for long long
  JSON::object::operator long long() const {
    switch (type) {
      {
        case object_types::INTEGER: {
          const std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          return *content;
        }
        case object_types::FLOAT: {
          const std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          if (*content > static_cast<double>(std::numeric_limits<long long>::max())) {
            return std::numeric_limits<long long>::max();
          }
          if (*content < static_cast<double>(std::numeric_limits<long long>::lowest())) {
            return std::numeric_limits<long long>::lowest();
          }
          return static_cast<long long>(*content);
        }
        case object_types::T: return 1;
        case object_types::F: return 0;
        default: throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON long long conversion failed: element is not a number: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief conversion function for unsigned long long
  JSON::object::operator unsigned long long() const {
    const long long result = static_cast<long long>(*this);
    if (result >= 0) {
      return result;
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON unsigned long long conversion failed: element is not a positive number: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
  }

  ///@brief conversion function for unsigned long
  JSON::object::operator unsigned long() const {
    return static_cast<unsigned long long>(*this);
  }

  ///@brief conversion function for int
  JSON::object::operator unsigned int() const {
    return static_cast<unsigned long long>(*this);
  }

  ///@brief conversion function for long
  JSON::object::operator long() const {
    return static_cast<long long>(*this);
  }

  ///@brief conversion function for int
  JSON::object::operator int() const {
    return static_cast<long long>(*this);
  }

  ///@brief conversion function for char
  JSON::object::operator char() const {
    return static_cast<long long>(*this);
  }

  ///@brief conversion function for string
  JSON::object::operator std::string() const {
    return this->toString(false, false);
  }

  ///@brief conversion function for vector of JSON::object
  JSON::object::operator std::vector<object>() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::vector conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    return *content;
  }

  ///@brief conversion function for map of JSON::object
  JSON::object::operator std::map<std::string, JSON::object>() const {
    if (type != object_types::DICTIONARY) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    return *content;
  }

  ///@brief allow to append to JSON::object_type::LIST
  void JSON::object::push_back(const JSON::object &content) const {
    if (type == JSON::object_types::LIST) {
      const std::shared_ptr<JSON::List> list_obj = std::static_pointer_cast<JSON::List>(obj);
      list_obj->push_back(content);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "push_back is just allowed for a JSON LIST", _VALUE_ILLEGAL_);
    }
  }

  ///@brief allow to merge two dictionaries or lists
  ///@note keys already stored in the dictionary will be updated
  void JSON::object::join(const JSON::object &content) const {
    if (type == JSON::object_types::DICTIONARY and content.type == JSON::object_types::DICTIONARY) {
      const std::shared_ptr<JSON::Dictionary> dict_obj = std::static_pointer_cast<JSON::Dictionary>(obj);
      const std::shared_ptr<JSON::Dictionary> dict_content = std::static_pointer_cast<JSON::Dictionary>(content.obj);
      for (auto &[key, value] : *dict_content) {
        dict_obj->operator[](key) = value;
      }
    } else if (type == JSON::object_types::LIST and content.type == JSON::object_types::LIST) {
      const std::shared_ptr<JSON::List> list_obj = std::static_pointer_cast<JSON::List>(obj);
      const std::shared_ptr<JSON::List> list_content = std::static_pointer_cast<JSON::List>(content.obj);
      for (auto &it : *list_content) {
        list_obj->push_back(it);
      }
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "insert is just allowed for a JSON DICTIONARY and LIST", _VALUE_ILLEGAL_);
    }
  }

  ///@brief gives the size of JSON::object_tpe::LIST, JSON::object_tpe::DICTIONARY or JSON::object_tpe::STRING
  size_t JSON::object::size() const {
    switch (type) {
      case object_types::DICTIONARY: {
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        return content->size();
      }
      case object_types::LIST: {
        const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
        return content->size();
      }
      case object_types::STRING: {
        const std::shared_ptr<std::string> content = std::static_pointer_cast<std::string>(obj);
        return content->size();
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "size is just allowed for a JSON DICTIONARY, LIST, or STRING", _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief checks if JSON::object_tpe::LIST, JSON::object_tpe::DICTIONARY or JSON::object_tpe::STRING is empty
  bool JSON::object::empty() const {
    if ((type == object_types::DICTIONARY) || (type == object_types::LIST) || (type == object_types::STRING)) {
      return not static_cast<bool>(*this);
    }
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "empty is just allowed for a JSON DICTIONARY, LIST, or STRING", _VALUE_ILLEGAL_);
  }

  ///@brief count the occurrence of key string in a JSON::object_type::DICTIONARY
  size_t JSON::object::count(const std::string &key) const {
    switch (type) {
      case object_types::DICTIONARY: {
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        return content->count(key);
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "count is just allowed for a JSON DICTIONARY", _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief returns an iterator to the end() of a JSON::object_type::DICTIONARY
  std::map<std::string, JSON::object>::iterator JSON::object::end() const {
    switch (type) {
      case object_types::DICTIONARY: {
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        return content->end();
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "end is just allowed for a JSON DICTIONARY", _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief returns an iterator to the begin() of a JSON::object_type::DICTIONARY
  std::map<std::string, JSON::object>::iterator JSON::object::begin() const {
    switch (type) {
      case object_types::DICTIONARY: {
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        return content->begin();
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "begin is just allowed for a JSON DICTIONARY", _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief returns an iterator to the result of find in JSON::object_type::DICTIONARY
  std::map<std::string, JSON::object>::iterator JSON::object::find(const std::string &key) const {
    switch (type) {
      case object_types::DICTIONARY: {
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        return content->find(key);
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "find is just allowed for a JSON DICTIONARY", _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief returns a vector of all keys in a JSON::object_type::DICTIONARY
  std::vector<std::string> JSON::object::keys() const {
    switch (type) {
      case object_types::DICTIONARY: {
        const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        std::vector<std::string> keys;
        keys.reserve(content->size());
        for (const auto &[key, value] : *content) {
          keys.push_back(key);
        }
        return keys;
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "keys is just allowed for a JSON DICTIONARY", _VALUE_ILLEGAL_);
      }
    }
  }

} // namespace aurostd

namespace aurostd {
  /// @namespace aurostd::JSON
  /// @brief unified namespace to read and write JSON
  ///
  /// @authors
  /// @mod{AS,2020,created JSONWriter}
  /// @mod{HE,20220924,rewrite to include parsing}
  ///
  /// @see
  /// @xlink{"JSON definition",https://www.json.org/json-en.html}
  /// @xlink{"Parsing JSON is a Minefield",https://seriot.ch/projects/parsing_json.html}
  ///
  /// Load data example
  /// @code
  /// aurostd::JSON::object jo;
  /// jo = aurostd::JSON::loadFile("testcases.json");
  /// // output complete file as JSON
  /// cout << jo.toString() << endl;
  /// // output element (alternative to .toSring())
  /// cout << (string) jo["xvector"][3] << endl;
  /// // save element
  /// xvector<double> xvd = jo["xvector"];
  /// cout << xvd << endl;
  /// // save sub element
  /// uint el_uint = jo["xvector"][3];
  /// cout << el_uint << endl;
  /// // cast to type
  /// cout << (float) jo["xvector"][3] << endl;
  /// // iterate over Dictionary
  /// for (const auto & [key, value]: static_cast<aurostd::JSON::Dictionary>(jo["a_dict"])) {
  ///   cout << key << " | " << (string) value << endl;
  /// }
  /// // iterate over List
  /// for (const auto & value: static_cast<aurostd::JSON::List>(jo["a_list"])) {
  ///   cout << (string) value << endl;
  /// }
  ///
  /// @endcode
  ///
  /// Save data example
  /// @code
  /// aurostd::JSON::object jo(aurostd::JSON::object_types::DICTIONARY);
  /// jo["string"] = "Hello World";
  /// jo["number"] = 4562318;
  /// jo["null"] = nullptr;
  /// jo["unicode"] = "🎃";
  /// jo["vector"] = (vector<float>){2.3,5.6,34.6};
  /// jo["xvector"] = (vector<double>){9.634,4.6,1E20};
  /// cout << jo.toString() << endl;
  /// jo.saveFile("testme.json");
  /// // or
  /// aurostd::JSON::saveFile(jo, "testme.json");
  /// @endcode

  /// @brief find a char inbetween a boundary
  /// @param content_ptr string pointer (string.c_str())
  /// @param border search boundary
  /// @param to_find char to find
  /// @return position
  static inline size_t range_find(const char *content_ptr, const std::pair<size_t, size_t> &border, const char to_find) {
    return static_cast<const char *>(memchr(content_ptr + border.first, to_find, border.second - border.first + 1)) - content_ptr;
  }

  /// @brief parse JSON string
  /// @param raw_content full json string
  /// @param border string borders
  /// @return c++ string
  /// @authors
  /// @mod{HE,20220924,created}
  std::string JSON::parse_string(const std::string &raw_content, std::pair<size_t, size_t> border) {
    if (border.second == 0) {
      border.second = raw_content.size() - 1;
    }
    size_t last_escape_pos = border.first;
    size_t current_escape_pos = range_find(raw_content.c_str(), border, '\\');

    if (current_escape_pos > border.second) {
      return raw_content.substr(border.first, border.second - border.first + 1);
    }

    std::string result;
    while (current_escape_pos <= border.second) {
      result += raw_content.substr(last_escape_pos, current_escape_pos - last_escape_pos);
      switch (raw_content[current_escape_pos + 1]) {
        case '"': {
          result += '"';
        } break;
        case '\\': {
          result += '\\';
        } break;
        case '/': {
          result += '/';
        } break;
        case 'b': {
          result += '\b';
        } break;
        case 'f': {
          result += '\f';
        } break;
        case 'n': {
          result += '\n';
        } break;
        case 'r': {
          result += '\r';
        } break;
        case 't': {
          result += '\t';
        } break;
        case 'u': {
          if (current_escape_pos + 6 > raw_content.size()) {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
          }
          result += unescape_unicode(raw_content, current_escape_pos);
          current_escape_pos += 2;
        } break;
        default: {
          stringstream message;
          message << "JSON parsing failed: string contains undefined escape '\\" << raw_content[current_escape_pos + 1] << "'";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message.str(), _FILE_WRONG_FORMAT_);
        }
      }
      last_escape_pos = current_escape_pos + 2;
      current_escape_pos = range_find(raw_content.c_str(), {last_escape_pos, border.second}, '\\');
    }
    result += raw_content.substr(last_escape_pos, border.second - last_escape_pos + 1);
    return result;
  }

  ///@brief escape characters to JSON
  std::string JSON::char_escape(const char16_t c) {
    switch (c) {
      default: return "";
      case ('"'): return "\\\"";
      case ('\\'): return "\\\\";
      case ('/'): return "\\/";
      case ('\b'): return "\\b";
      case ('\f'): return "\\f";
      case ('\n'): return "\\n";
      case ('\r'): return "\\r";
      case ('\t'): return "\\t";
    }
  }

  /// @brief prepare string for JSON output with or without Unicode escapes
  /// @param raw string to escape
  /// @param unicode if true escape Unicode characters (default true)
  /// @return JSON string
  /// @authors
  /// @mod{HE,20221109,created}
  std::string JSON::escape(const std::string &raw, const bool unicode) {
    std::stringstream out;
    if (unicode) {
      const std::u16string utf16 = std::wstring_convert<std::codecvt_utf8_utf16<char16_t>, char16_t>{}.from_bytes(raw.data());
      for (const char16_t c : utf16) {
        if (c < 0x80) {
          if (char_escape(c).empty()) {
            out << static_cast<char>(c);
          } else {
            out << char_escape(c);
          }
        } else {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << c;
        }
      }
    } else {
      for (const char16_t c : raw) {
        if (char_escape(c).empty()) {
          out << static_cast<char>(c);
        } else {
          out << char_escape(c);
        }
      }
    }
    return out.str();
  }

  /// @brief convert a unicode codepoint to a series of utf8 chars
  /// @param cp 32bit codepoint
  /// @return utf8 representation
  /// @authors
  /// @mod{HE,20221109,created}
  std::string JSON::char32_to_string(const char32_t cp) {
    string out;
    if (cp < 0x80) {
      out += static_cast<char>(cp);
      return out;
    }
    if (cp < 0x800) {
      out += static_cast<char>((cp >> 6) | 0xc0);
      out += static_cast<char>((cp & 0x3f) | 0x80);
      return out;
    }
    if (cp < 0x10000) {
      out += static_cast<char>((cp >> 12) | 0xe0);
      out += static_cast<char>(((cp >> 6) & 0x3f) | 0x80);
      out += static_cast<char>((cp & 0x3f) | 0x80);
      return out;
    }

    out += static_cast<char>((cp >> 18) | 0xf0);
    out += static_cast<char>(((cp >> 12) & 0x3f) | 0x80);
    out += static_cast<char>(((cp >> 6) & 0x3f) | 0x80);
    out += static_cast<char>((cp & 0x3f) | 0x80);
    return out;
  }

  /// @brief unescape JSON unicode instances
  /// @param raw raw content string
  /// @param pos starting position of the hex representation
  /// @return unescaped utf8 characters
  /// @note this function supports pairs
  /// @authors
  /// @mod{HE,20221109,created}
  std::string JSON::unescape_unicode(const std::string &raw, size_t &pos) {
    pos += 2;
    const char32_t cp = (char32_t) aurostd::string2utype<uint>(raw.substr(pos, 4), 16);

    if (cp < 0xd800 || cp > 0xdfff) {
      return char32_to_string(cp);
    } else if (cp > 0xdbff) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
    } else {
      if (raw[pos + 4] == '\\' and raw[pos + 5] == 'u') {
        pos += 6;
        const char32_t trailing_cp = static_cast<char32_t>(aurostd::string2utype<uint>(raw.substr(pos, 4), 16));
        if (trailing_cp < 0xdc00 || trailing_cp > 0xdfff) {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
        }
        const char32_t combo_cp = ((cp - 0xd800) << 10) + (trailing_cp - 0xdc00) + 0x10000;
        return char32_to_string(combo_cp);
      } else {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
      }
    }
  }

  /// @brief strip whitespaces
  /// @param raw_content full json string
  /// @param border working border
  /// @return stripped raw string borders
  /// @authors
  /// @mod{HE,20220924,created}
  std::pair<size_t, size_t> JSON::find_strip(const std::string &raw_content, std::pair<size_t, size_t> border) {
    if (border.second == 0) {
      border.second = raw_content.size() - 1;
    }
    const size_t start = raw_content.find_first_not_of(" \n\t\r\v\f", border.first);
    const size_t end = raw_content.find_last_not_of(" \n\t\r\v\f", border.second);
    return {start, min(end, border.second)};
  }

  /// @brief find the border of a JSON string
  /// @param raw_content full json string
  /// @param border working border
  /// @return string borders
  /// @authors
  /// @mod{HE,20220924,created}
  std::pair<size_t, size_t> JSON::find_string(const std::string &raw_content, std::pair<size_t, size_t> border) {
    if (border.second == 0) {
      border.second = raw_content.size() - 1;
    }
    const size_t start = range_find(raw_content.c_str(), border, '"');
    if (start >= border.second) {
      return {std::string::npos, std::string::npos};
    }
    size_t end = start;
    uint escape_check = 0;
    size_t escape_check_pos = 0;
    if (start == std::string::npos) {
      return {start, end};
    }
    do {
      end = range_find(raw_content.c_str(), {end + 1, border.second}, '"');
      escape_check_pos = end - 1;
      escape_check = 0;
      while (raw_content[escape_check_pos] == '\\' and escape_check_pos >= start) {
        escape_check++;
        escape_check_pos--;
      }
    } while ((escape_check % 2 != 0) && (end <= border.second));
    if (end > border.second) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: string not closed", _FILE_WRONG_FORMAT_);
    }
    return {start, end};
  }

  /// @brief find the border of an encapsulated by brackets
  /// @param raw_content full json string
  /// @param kind_open type of bracket [ { (
  /// @param border working border
  /// @return bracket borders
  /// @authors
  /// @mod{HE,20220924,created}
  std::pair<size_t, size_t> JSON::find_bracket(const std::string &raw_content, const char kind_open, std::pair<size_t, size_t> border) {
    char kind_close;
    switch (kind_open) {
      case '[': kind_close = ']'; break;
      case '{': kind_close = '}'; break;
      case '(': kind_close = ')'; break;
      default: throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: bracket kind not allowed", _FILE_WRONG_FORMAT_);
    }
    const size_t start = raw_content.find(kind_open, border.first);
    size_t end = start;
    if (start > border.second) {
      return {std::string::npos, std::string::npos};
    }
    size_t next_open = 0;
    size_t next_close = 0;
    std::pair<size_t, size_t> string_section;
    size_t open_count = 1;
    do {
      // cerr <<"end: " << end << " | "  << "open: " << open_count << " | current selection: " << raw_content.substr(start, end-start+1) << endl; // detailed debug
      string_section = find_string(raw_content, {end + 1, border.second});
      next_close = range_find(raw_content.c_str(), {end + 1, border.second}, kind_close);
      next_open = range_find(raw_content.c_str(), {end + 1, border.second}, kind_open);
      if (next_close > string_section.first && next_open > string_section.first) {
        end = string_section.second;
        continue;
      }
      if (next_close < next_open) {
        if (next_close < string_section.first) {
          end = next_close;
          open_count--;
        }
      } else if (next_close > next_open) {
        if (next_open < string_section.first) {
          end = next_open;
          open_count++;
        }
      } else { // when == std::string::npos
        break;
      }

    } while (open_count > 0 && (end < border.second));

    return {start, end};
  }

  /// @brief parse a raw JSON string
  /// @param raw_content full json string
  /// @param border working border
  /// @return parsed json object
  /// @authors
  /// @mod{HE,20220924,created}
  JSON::object JSON::parse(const std::string &raw_content, std::pair<size_t, size_t> border) {
    if (border.second == 0) {
      border.second = raw_content.size() - 1;
    }
    object result;
    result.type = object_types::NONE;
    result.obj = nullptr;
    border = find_strip(raw_content, border);
    std::pair<size_t, size_t> section({0, 0});
    switch (raw_content[border.first]) {
      case '{': {
        if (raw_content[border.second] != '}') {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: object not closed by '}'", _FILE_WRONG_FORMAT_);
        }
        const std::shared_ptr<JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
        border = {border.first + 1, border.second - 1};
        section = {0, 0};
        do {
          border = find_strip(raw_content, border);
          // cout  << border.first << " | "<< border.second << " - current: " << raw_content.substr(border.first, border.second-border.first+1) << endl; // detailed debug
          if (border.first >= border.second) {
            break;
          }
          section = find_string(raw_content, border);
          if (section.first == std::string::npos) {
            break;
          }
          const std::string key_string = parse_string(raw_content, {section.first + 1, section.second - 1}); // cut " already
          //  cout << "key: " << key_string << endl;
          border.first = raw_content.find(':', section.second) + 1;
          border = find_strip(raw_content, border);
          // cout  << border.first << " | "<< border.second << " - switch: " << raw_content.substr(border.first, border.second-border.first+1) << endl; // detailed debug

          switch (raw_content[border.first]) {
            case '"': {
              section = find_string(raw_content, border);
            } break;
            case '[': {
              section = find_bracket(raw_content, '[', border);
            } break;
            case '{': {
              section = find_bracket(raw_content, '{', border);
            } break;
            default: {
              section.first = border.first;
              section.second = min(border.second, range_find(raw_content.c_str(), border, ',') - 1);
            }
          }
          new_dictionary->insert({key_string, parse(raw_content, section)});
          border.first = range_find(raw_content.c_str(), {section.second, border.second}, ',');
          if (border.first == std::string::npos) {
            border.first = border.second;
          } else {
            border.first++;
          }
        } while (border.first < border.second);
        result.obj = new_dictionary;
        result.type = object_types::DICTIONARY;
      } break;
      case '[': {
        if (raw_content[border.second] != ']') {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: list not closed by ']'", _FILE_WRONG_FORMAT_);
        }
        const std::shared_ptr<JSON::List> new_list = std::make_shared<JSON::List>();
        border = {border.first + 1, border.second - 1};
        section = {0, 0};
        do {
          border = find_strip(raw_content, border);
          switch (raw_content[border.first]) {
            case '"': {
              section = find_string(raw_content, border);
            } break;
            case '[': {
              section = find_bracket(raw_content, '[', border);
            } break;
            case '{': {
              section = find_bracket(raw_content, '{', border);
            } break;
            default: {
              section.first = border.first;
              section.second = min(border.second, raw_content.find(',', border.first) - 1);
            }
          }
          new_list->emplace_back(parse(raw_content, section));
          border.first = raw_content.find(',', section.second);
          if (border.first == std::string::npos) {
            border.first = border.second;
          } else {
            border.first++;
          }
        } while (section.second < border.second);
        result.obj = new_list;
        result.type = object_types::LIST;
      } break;
      case '"': {
        if (raw_content[border.second] != '"') {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: string not enclosed by '\"'", _FILE_WRONG_FORMAT_);
        }
        const std::shared_ptr<std::string> new_string = std::make_shared<std::string>();
        *new_string = parse_string(raw_content, {border.first + 1, border.second - 1});
        result.obj = new_string;
        result.type = object_types::STRING;
      } break;
      case 'n': {
        result.obj = nullptr;
        result.type = object_types::NONE;
      } break;
      case 't': {
        result.obj = nullptr;
        result.type = object_types::T;
      } break;
      case 'f': {
        result.obj = nullptr;
        result.type = object_types::F;
      } break;
      default: { // number
        if (raw_content.find('.', border.first) != std::string::npos || raw_content.find('E', border.first) != std::string::npos || raw_content.find('e', border.first) != std::string::npos) {
          const std::shared_ptr<double> new_number = std::make_shared<double>();
          *new_number = std::strtod(raw_content.c_str() + border.first, nullptr);
          result.obj = new_number;
          result.type = object_types::FLOAT;
        } else {
          const std::shared_ptr<long long> new_number = std::make_shared<long long>();
          *new_number = std::strtoll(raw_content.c_str() + border.first, nullptr, 10);
          result.obj = new_number;
          result.type = object_types::INTEGER;
        }
      }
    }
    return result;
  }

  /// @brief create a JSON::object from file
  /// @param file_path file path to load from
  /// @return parsed JSON::object
  JSON::object JSON::loadFile(const std::string &file_path) {
    const std::string raw_content = aurostd::file2string(file_path);
    return parse(raw_content);
  }

  /// @brief create a JSON::object from raw string
  /// @param content raw JSON content
  /// @return parsed JSON::object
  JSON::object JSON::loadString(const std::string &content) {
    return parse(content);
  }

  /// @brief convert JSON::object to string
  /// @param root JSON::object to convert
  /// @param escape_unicode if true escape Unicode characters
  /// @return converted string
  std::string JSON::toString(const object &root, const bool escape_unicode) {
    return root.toString(true, escape_unicode);
  }

  /// @brief save JSON::object to file
  /// @param root JSON::object to save
  /// @param file_path file path to save to
  /// @param escape_unicode if true escape Unicode characters
  void JSON::saveFile(const object &root, const std::string &file_path, const compression_type ct, const bool escape_unicode) {
    aurostd::string2file(root.toString(true, escape_unicode), file_path, ct);
  }

} // namespace aurostd

#endif // _AUROSTD_XPARSER_JSON_CPP_
