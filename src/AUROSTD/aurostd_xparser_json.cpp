// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// This JSON class is the evolution of different prior solutions to integrate JSON with the AFLOW source base.
// hagen.eckert@duke.edu


#ifndef _AUROSTD_XPARSER_JSON_CPP_
#define _AUROSTD_XPARSER_JSON_CPP_

#ifndef _AUROSTD_XPARSER_H_

#include "aurostd_xparser.h"

#endif

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
  JSON::object &JSON::object::operator[](const size_t index) const{
    if (this->type == JSON::object_types::LIST) {
      std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(this->obj);
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
  JSON::object &JSON::object::operator[](const std::string key) const{
    if (this->type == JSON::object_types::DICTIONARY) {
      std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(this->obj);
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
      std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(this->obj);
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
        std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        for (const auto &entry: *content) {
          if (first) first = false;
          else result << ",";
          result << "\"" << entry.first << "\":" << entry.second; //
        }
        result << "}";
      }
        break;
      case object_types::LIST: {
        result << "[";
        std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
        for (const auto &entry: *content) {
          if (first) first = false;
          else result << ",";
          result << entry; //
        }
        result << "]";
      }
        break;
      case object_types::STRING: {
        std::shared_ptr <std::string> content = std::static_pointer_cast<std::string>(obj);
        if (json_format) result << "\"" << JSON::escape(*content, escape_unicode) << "\"";
        else return *content;
      }
        break;
      case object_types::INTEGER: {
        std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
        result << *content;
      }
        break;
      case object_types::FLOAT: {
        std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
        result << *content;
      }
        break;
      case object_types::T: {
        result << "true";
      }
        break;
      case object_types::F: {
        result << "false";
      }
        break;
      case object_types::NONE: {
        result << "null";
      }
        break;
      default:
        break;
    }
    return result.str();
  }

  /// @brief save JSON::object to file
  /// @param file_path path to save to
  /// @param escape_unicode if `true` escape unicode in strings - `αβγ` becomes `\u03b1\u03b2\u03b3` (default: `true`)
  /// @authors
  /// @mod{HE,20221110,created}
  void JSON::object::saveFile(const std::string & file_path, const bool escape_unicode) const {
    aurostd::string2file(this->toString(true, escape_unicode), file_path);
  }

  /// @brief change this JSON::object to a JSON::object_types::STRING
  /// @authors
  /// @mod{HE,20220926,created}
  void JSON::object::fromString(const std::string & content){
    std::shared_ptr <std::string> new_string = std::make_shared<std::string>();
    *new_string = content;
    this->obj = new_string;
    this->type = object_types::STRING;
  }

  /// @brief change this JSON::object to a JSON::object_types::INTEGER
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype>  void JSON::object::fromNumber(const utype content){
    if (std::is_integral<utype>::value){
      this->type = object_types::INTEGER;
      std::shared_ptr <long long> new_content = std::make_shared<long long>();
      *new_content = (long long) content;
      this->obj = new_content;
      return;
    } else if (std::is_floating_point<utype>::value){
      if (std::isnan(content)) {
        this->type = object_types::NONE;
        this->obj = nullptr;
        return;
      }
      this->type = object_types::FLOAT;
      std::shared_ptr <double> new_content = std::make_shared<double>();
      *new_content = (double) content;
      this->obj = new_content;
      return;
    }
  }

  /// @brief change this JSON::object to a complex representation
  /// @authors
  /// @mod{HE,20221118,created}
  template<class utype>  void JSON::object::fromComplex(const xcomplex<utype> content){
    std::shared_ptr <JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
    this->obj = new_dictionary;
    this->type = object_types::DICTIONARY;
    new_dictionary->insert({"real", content.re});
    new_dictionary->insert({"imag", content.im});
    return;
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a xvector
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromXvector(const xvector<utype> &content) {
    std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (int idx=content.lrows; idx<=content.urows; idx++) {
      new_list->push_back(content[idx]);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a xmatrix
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromXmatrix(const xmatrix<utype> &content) {
    std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (int idx=content.lrows; idx<=content.urows; idx++) {
      new_list->push_back(content(idx));
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a vector
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromVector(const vector<utype> & content) {
    std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (utype entry: content) {
      new_list->push_back(entry);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::DICTIONARY based on a map
  /// @authors
  /// @mod{HE,20220926,created}
  template<class utype> void JSON::object::fromMap(const std::map<std::string, utype> & content) {
    std::shared_ptr <JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
    this->obj = new_dictionary;
    this->type = object_types::DICTIONARY;
    for (std::pair<std::string, utype> entry: content) {
      new_dictionary->insert({entry.first, entry.second});
    }
  }

  /// @brief converting constructor: set the content of this JSON::object based on a char
  JSON::object::object(const char* content) {
    this->fromString(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a string
  JSON::object::object(const std::string &content) {
    this->fromString(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a bool
  JSON::object::object(bool content) {
    if (content) this->type = object_types::T;
    else this->type = object_types::F;
  }

  /// @brief converting constructor: set the content of this JSON::object to null
  /// @note pass in a nullptr
  JSON::object::object(std::nullptr_t content) {
    this->type = object_types::NONE;
    this->obj = content;
  }

  /// @brief converting constructor: set the content of this JSON::object based on a number
  template<class utype> JSON::object::object(const utype content) {
    this->fromNumber(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a number
  template<class utype> JSON::object::object(const xcomplex<utype> &content) {
    this-> fromComplex(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a xvector
  template<class utype> JSON::object::object(const xvector<utype> & content) {
    this->fromXvector(content);
  }
  // template instantiation for xvector types
  template JSON::object::object(const xvector<char> &);
  template JSON::object::object(const xvector<int> &);
  template JSON::object::object(const xvector<unsigned int> &);
  template JSON::object::object(const xvector<long long> &);
  template JSON::object::object(const xvector<unsigned long long> &);
  template JSON::object::object(const xvector<float> &);
  template JSON::object::object(const xvector<double> &);
  template JSON::object::object(const xvector<xcomplex<float>> &);
  template JSON::object::object(const xvector<xcomplex<double>> &);

  /// @brief converting constructor: set the content of this JSON::object based on a xmatrix
  template<class utype> JSON::object::object(const xmatrix<utype> & content) {
    this->template fromXmatrix(content);
  }
  // template instantiation for xmatrix types
  template JSON::object::object(const xmatrix<int> &);
  template JSON::object::object(const xmatrix<unsigned int> &);
  template JSON::object::object(const xmatrix<long long> &);
  template JSON::object::object(const xmatrix<unsigned long long> &);
  template JSON::object::object(const xmatrix<float> &);
  template JSON::object::object(const xmatrix<double> &);
  template JSON::object::object(const xmatrix<xcomplex<float>> &);
  template JSON::object::object(const xmatrix<xcomplex<double>> &);

  /// @brief converting constructor: set the content of this JSON::object based on a vector
  template<class utype> JSON::object::object(const vector<utype> & content) {
    this->fromVector(content);
  }
  // template instantiation for vector types
  template JSON::object::object(const std::vector<int> &);
  template JSON::object::object(const std::vector<unsigned int> &);
  template JSON::object::object(const std::vector<long long> &);
  template JSON::object::object(const std::vector<unsigned long long> &);
  template JSON::object::object(const std::vector<float> &);
  template JSON::object::object(const std::vector<double> &);
  template JSON::object::object(const std::vector<xcomplex<float>> &);
  template JSON::object::object(const std::vector<xcomplex<double>> &);
  template JSON::object::object(const std::vector<std::string> &);
  template JSON::object::object(const std::vector<std::vector<int>> &);
  template JSON::object::object(const std::vector<std::vector<unsigned int>> &);
  template JSON::object::object(const std::vector<std::vector<long long>> &);
  template JSON::object::object(const std::vector<std::vector<unsigned long long>> &);
  template JSON::object::object(const std::vector<std::vector<float>> &);
  template JSON::object::object(const std::vector<std::vector<double>> &);
  template JSON::object::object(const std::vector<std::vector<xcomplex<float>>> &);
  template JSON::object::object(const std::vector<std::vector<xcomplex<double>>> &);
  template JSON::object::object(const std::vector<std::vector<std::string>> &);

  /// @brief converting constructor: set the content of this JSON::object based on a map
  template<class utype> JSON::object::object(const std::map<std::string, utype> & content) {
    this->template fromMap(content);
  }
  // template instantiation for map types (key is always a sting in JSON)
  template JSON::object::object(const std::map<std::string, int> &);
  template JSON::object::object(const std::map<std::string, unsigned int> &);
  template JSON::object::object(const std::map<std::string, long long> &);
  template JSON::object::object(const std::map<std::string, unsigned long long> &);
  template JSON::object::object(const std::map<std::string, float> &);
  template JSON::object::object(const std::map<std::string, double> &);
  template JSON::object::object(const std::map<std::string, xcomplex<float>> &);
  template JSON::object::object(const std::map<std::string, xcomplex<double>> &);
  template JSON::object::object(const std::map<std::string, std::string> &);
  template JSON::object::object(const std::map<std::string, std::vector<float>> &);

  /// @brief create an empty JSON::object of a given JSON::object_types
  /// @authors
  /// @mod{HE,20221031,created}
  JSON::object::object(object_types create_type) {
    this->type = create_type;
    switch (create_type) {
      {
        case object_types::DICTIONARY: {
          std::shared_ptr <JSON::Dictionary> content = std::make_shared<JSON::Dictionary>();
          this->obj = content;
          break;
        }
        case object_types::LIST: {
          std::shared_ptr <JSON::List> content = std::make_shared<JSON::List>();
          this->obj = content;
          break;
        }
        case object_types::STRING: {
          std::shared_ptr <std::string> content = std::make_shared<std::string>();
          this->obj = content;
          break;
        }
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::make_shared<long long int>();
          this->obj = content;
          break;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::make_shared<double>();
          this->obj = content;
          break;
        }
        default: {
          this->obj = nullptr;
        }
      }
    }
  }

  ///@brief assignment operator for char
  void JSON::object::operator=(const char* content) {
    this->fromString(content);
  }

  ///@brief assignment operator for string
  void JSON::object::operator=(const std::string & content) {
    this->fromString(content);
  }

  ///@brief assignment operator for bool
  void JSON::object::operator=(bool content) {
      if (content) this->type = object_types::T;
      else this->type = object_types::F;
  }

  ///@brief assignment operator for nullptr
  void JSON::object::operator=(std::nullptr_t content) {
    this->type = object_types::NONE;
    this->obj = content;
  }

  ///@brief assignment operator for numbers
  template<class utype> void JSON::object::operator=(const utype content){
    this->template fromNumber(content);
  }
  // template instantiation for number types
  template void JSON::object::operator=(const char);
  template void JSON::object::operator=(const int);
  template void JSON::object::operator=(const unsigned int);
  template void JSON::object::operator=(const long long);
  template void JSON::object::operator=(const unsigned long long);
  template void JSON::object::operator=(const double);
  template void JSON::object::operator=(const float);

  ///@brief assignment operator for xcomplex
  template<class utype> void JSON::object::operator=(const xcomplex<utype> & content){
    fromComplex(content);
  }
  template void JSON::object::operator=(const xcomplex<float> &);
  template void JSON::object::operator=(const xcomplex<double> &);

  ///@brief assignment operator for xvector
  template<class utype> void JSON::object::operator=(const xvector<utype> & content){
    fromXvector(content);
  }
  // template instantiation for xvector types
  template void JSON::object::operator=(const xvector<int> &);
  template void JSON::object::operator=(const xvector<unsigned int> &);
  template void JSON::object::operator=(const xvector<long long> &);
  template void JSON::object::operator=(const xvector<unsigned long long> &);
  template void JSON::object::operator=(const xvector<float> &);
  template void JSON::object::operator=(const xvector<double> &);
  template void JSON::object::operator=(const xvector<xcomplex<float>> &);
  template void JSON::object::operator=(const xvector<xcomplex<double>> &);

  ///@brief assignment operator for xmatrix
  template<class utype> void JSON::object::operator=(const xmatrix<utype> & content){
    fromXmatrix(content);
  }
  // template instantiation for xmatrix types
  template void JSON::object::operator=(const xmatrix<int> &);
  template void JSON::object::operator=(const xmatrix<unsigned int> &);
  template void JSON::object::operator=(const xmatrix<long long> &);
  template void JSON::object::operator=(const xmatrix<unsigned long long> &);
  template void JSON::object::operator=(const xmatrix<float> &);
  template void JSON::object::operator=(const xmatrix<double> &);
  template void JSON::object::operator=(const xmatrix<xcomplex<float>> &);
  template void JSON::object::operator=(const xmatrix<xcomplex<double>> &);

  ///@brief assignment operator for vector
  template<class utype> void JSON::object::operator=(const std::vector<utype> & content){
    fromVector(content);
  }
  // template instantiation for vector types
  template void JSON::object::operator=(const std::vector<int> &);
  template void JSON::object::operator=(const std::vector<unsigned int> &);
  template void JSON::object::operator=(const std::vector<long long> &);
  template void JSON::object::operator=(const std::vector<unsigned long long> &);
  template void JSON::object::operator=(const std::vector<float> &);
  template void JSON::object::operator=(const std::vector<double> &);
  template void JSON::object::operator=(const std::vector<xcomplex<float>> &);
  template void JSON::object::operator=(const std::vector<xcomplex<double>> &);
  template void JSON::object::operator=(const std::vector<std::string> &);
  template void JSON::object::operator=(const std::vector<vector<int>> &);
  template void JSON::object::operator=(const std::vector<vector<unsigned int>> &);
  template void JSON::object::operator=(const std::vector<vector<long long>> &);
  template void JSON::object::operator=(const std::vector<vector<unsigned long long>> &);
  template void JSON::object::operator=(const std::vector<vector<float>> &);
  template void JSON::object::operator=(const std::vector<vector<double>> &);
  template void JSON::object::operator=(const std::vector<vector<xcomplex<float>>> &);
  template void JSON::object::operator=(const std::vector<vector<xcomplex<double>>> &);
  template void JSON::object::operator=(const std::vector<vector<std::string>> &);
  template void JSON::object::operator=(const std::vector<xvector<int>> &);
  template void JSON::object::operator=(const std::vector<xvector<unsigned int>> &);
  template void JSON::object::operator=(const std::vector<xvector<long long>> &);
  template void JSON::object::operator=(const std::vector<xvector<unsigned long long>> &);
  template void JSON::object::operator=(const std::vector<xvector<float>> &);
  template void JSON::object::operator=(const std::vector<xvector<double>> &);

  ///@brief assignment operator for map
  template<class utype> void JSON::object::operator=(const std::map<std::string, utype> & content){
    fromMap(content);
  }
  // template instantiation for map types (key is always a sting in JSON)
  template void JSON::object::operator=(const std::map<std::string, int> &);
  template void JSON::object::operator=(const std::map<std::string, unsigned int> &);
  template void JSON::object::operator=(const std::map<std::string, long long> &);
  template void JSON::object::operator=(const std::map<std::string, unsigned long long> &);
  template void JSON::object::operator=(const std::map<std::string, float> &);
  template void JSON::object::operator=(const std::map<std::string, double> &);
  template void JSON::object::operator=(const std::map<std::string, xcomplex<float>> &);
  template void JSON::object::operator=(const std::map<std::string, xcomplex<double>> &);
  template void JSON::object::operator=(const std::map<std::string, std::string> &);
  template void JSON::object::operator=(const std::map<std::string, std::vector<float>> &);


  ///@brief conversion function for bool
  ///@note for LIST, DICTIONARY and STRING tests emptiness
  ///@note for NULL returns always false
  JSON::object::operator bool() const {
    switch (type) {
      {
        case object_types::DICTIONARY: {
          std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::LIST: {
          std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::STRING: {
          std::shared_ptr <std::string> content = std::static_pointer_cast<std::string>(obj);
          if (content->empty()) return false;
          else return true;
        }
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          if (*content) return true;
          else return false;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          if (*content) return true;
          else return false;
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
          std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          return (double) *content;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          return *content;
        }
        case object_types::T:
          return 1.0;
        case object_types::F:
          return 0.0;
        case object_types::NONE:
          return NAN;
        default:
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON double conversion failed: is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief conversion function for float
  JSON::object::operator float() const {
    return (double) *this;
  }

  ///@brief conversion function for long long
  JSON::object::operator long long() const {
    switch (type) {
      {
        case object_types::INTEGER: {
          std::shared_ptr<long long int> content = std::static_pointer_cast<long long int>(obj);
          return *content;
        }
        case object_types::FLOAT: {
          std::shared_ptr<double> content = std::static_pointer_cast<double>(obj);
          if (*content > (double) std::numeric_limits<long long>::max())
            return std::numeric_limits<long long>::max();
          if (*content < (double) std::numeric_limits<long long>::lowest())
            return std::numeric_limits<long long>::lowest();
          return (long long) *content;
        }
        case object_types::T:
          return 1;
        case object_types::F:
          return 0;
        default:
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON long long conversion failed: element is not a number: " + (std::string) * this, _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief conversion function for unsigned long long
  JSON::object::operator unsigned long long() const {
    long long result = (long long) *this;
    if (result >= 0) return result;
    else throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON unsigned long long conversion failed: element is not a positive number: " + (std::string) * this, _VALUE_ILLEGAL_);
  }

  ///@brief conversion function for unsigned long
  JSON::object::operator unsigned long() const {
    return (unsigned long long) *this;
  }

  ///@brief conversion function for int
  JSON::object::operator unsigned int() const {
    return (unsigned long long) *this;
  }

  ///@brief conversion function for long
  JSON::object::operator long() const {
    return (long long) *this;
  }

  ///@brief conversion function for int
  JSON::object::operator int() const {
    return (long long) *this;
  }

  ///@brief conversion function for int
  template<class utype> JSON::object::operator xcomplex<utype>() const {
    if (type != object_types::DICTIONARY)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xcomplex conversion failed: element is not a DICTIONARY: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    return xcomplex<utype>((utype) content->operator[]("real"), (utype) content->operator[]("imag"));
  }
  template JSON::object::operator xcomplex<double>() const;
  template JSON::object::operator xcomplex<float>() const;


  ///@brief conversion function for string
  JSON::object::operator std::string() const {
    return this->toString(false, false);
  }

  ///@brief conversion function for vector
  template<class utype> JSON::object::operator std::vector<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a LIST: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::vector <utype> result;
    std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) return result;
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (const JSON::object &so: *content)
        result.push_back((utype) so); // should be emplace_back - change in cpp14 (error for bool)
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (const JSON::object &so: *content) {
        if (so.type == object_types::DICTIONARY) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + (std::string) so, _VALUE_ILLEGAL_);
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + (std::string) so, _VALUE_ILLEGAL_);
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) result.push_back((utype) so);
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + (std::string) so, _VALUE_ILLEGAL_);
      }
    }
    return result;
  }
  // template instantiation for vector types
  template JSON::object::operator std::vector<std::string>() const;
  template JSON::object::operator std::vector<double>() const;
  template JSON::object::operator std::vector<float>() const;
  template JSON::object::operator std::vector<long long>() const;
  template JSON::object::operator std::vector<unsigned long long>() const;
  template JSON::object::operator std::vector<long>() const;
  template JSON::object::operator std::vector<int>() const;
  template JSON::object::operator std::vector<uint>() const;
  template JSON::object::operator std::vector<bool>() const;
  template JSON::object::operator std::vector<xcomplex<float>>() const;
  template JSON::object::operator std::vector<xcomplex<double>>() const;

  ///@brief conversion function for map of JSON::object
  JSON::object::operator std::map<std::string, JSON::object> () const {
    if (type != object_types::DICTIONARY)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    return *content;
  }

  ///@brief conversion function for map with type conversion
  template<class utype> JSON::object::operator std::map<std::string, utype>() const {
    if (type != object_types::DICTIONARY)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::map <std::string, utype> result;
    std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    if (content->empty()) return result;
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (std::pair<const std::string, object> &entry: *content)
        result.insert({entry.first, (utype) entry.second});
    } else if (std::is_floating_point<utype>::value) {
      for (std::pair<const std::string, object> &entry: *content) {
        if (entry.second.type > object_types::STRING) result.insert({entry.first, (utype) entry.second});
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + (std::string) entry.second, _VALUE_ILLEGAL_);
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (std::pair<const std::string, object> &entry: *content) {
        if (entry.second.type > object_types::STRING && entry.second.type < object_types::NONE)
          result.insert({entry.first, (utype) entry.second});
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + (std::string) entry.second, _VALUE_ILLEGAL_);
      }
    }
    return result;
  }
  // template instantiation for map types (key is always a sting in JSON)
  template JSON::object::operator std::map<std::string, std::string>() const;
  template JSON::object::operator std::map<std::string, double>() const;
  template JSON::object::operator std::map<std::string, float>() const;
  template JSON::object::operator std::map<std::string, long long>() const;
  template JSON::object::operator std::map<std::string, long>() const;
  template JSON::object::operator std::map<std::string, int>() const;
  template JSON::object::operator std::map<std::string, uint>() const;
  template JSON::object::operator std::map<std::string, bool>() const;

  ///@brief conversion function for xvector
  template<class utype> JSON::object::operator aurostd::xvector<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xvector conversion failed: element is not a LIST: " + (std::string) * this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    aurostd::xvector<utype> result(content->size(), 1);
    if (content->empty()) return result;
    size_t idx = 1;
    if (typeid(utype) == typeid(bool)) {
      for (const JSON::object &so: *content) {
        result[idx] = (utype) so;
        idx++;
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (const JSON::object &so: *content) {
        if (so.type == object_types::DICTIONARY) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + (std::string) so, _VALUE_ILLEGAL_);
        idx++;
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + (std::string) so, _VALUE_ILLEGAL_);
        idx++;
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so: *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) result[idx] = (utype) so;
        else
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + (std::string) so, _VALUE_ILLEGAL_);
        idx++;
      }
    }
    return result;
  }
  // template instantiation for xvector types
  template JSON::object::operator aurostd::xvector<double>() const;
  template JSON::object::operator aurostd::xvector<float>() const;
  template JSON::object::operator aurostd::xvector<long long>() const;
  template JSON::object::operator aurostd::xvector<unsigned long long>() const;
  template JSON::object::operator aurostd::xvector<long>() const;
  template JSON::object::operator aurostd::xvector<int>() const;
  template JSON::object::operator aurostd::xvector<uint>() const;
  template JSON::object::operator aurostd::xvector<bool>() const;
  template JSON::object::operator aurostd::xvector<xcomplex<double>>() const;
  template JSON::object::operator aurostd::xvector<xcomplex<float>>() const;

  ///@brief conversion function for xmatrix
  template<class utype> JSON::object::operator aurostd::xmatrix<utype>() const {
    if (type != object_types::LIST)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: element is not a LIST: " + (std::string) *this, _VALUE_ILLEGAL_);
    std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) return aurostd::xmatrix<utype>();

    // size scan
    size_t rows = 0;
    size_t cols = 0;
    std::vector <size_t> cols_v;
    for (const JSON::object &so: *content) {
      if (so.type != object_types::LIST)
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub element is not a LIST: " + (std::string) *this, _VALUE_ILLEGAL_);
      std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(so.obj);
      rows++;
      cols_v.template emplace_back(row->size());
    }
    if (std::adjacent_find(cols_v.begin(), cols_v.end(), std::not_equal_to<size_t>()) == cols_v.end()) cols = cols_v[0];
    else
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub elements have different length", _VALUE_ILLEGAL_);

    aurostd::xmatrix<utype> result(rows, cols);
    if (typeid(utype) == typeid(bool)) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          result[r][c] = (utype) row->operator[](c - 1);
        }
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type == object_types::DICTIONARY) result[r][c] = (utype) row->operator[](c - 1);
          else
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + (std::string) row->operator[](c - 1), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING) result[r][c] = (utype) row->operator[](c - 1);
          else
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + (string) row->operator[](c - 1), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (int r = result.lrows; r <= result.urows; r++) {
        std::shared_ptr <JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING && row->operator[](c - 1).type < object_types::NONE)
            result[r][c] = (utype) row->operator[](c - 1);
          else
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + (string) row->operator[](c - 1), _VALUE_ILLEGAL_);
        }
      }
    }
    return result;
  }
  // template instantiation for xmatrix types
  template JSON::object::operator aurostd::xmatrix<double>() const;
  template JSON::object::operator aurostd::xmatrix<float>() const;
  template JSON::object::operator aurostd::xmatrix<long long>() const;
  template JSON::object::operator aurostd::xmatrix<unsigned long long>() const;
  template JSON::object::operator aurostd::xmatrix<long>() const;
  template JSON::object::operator aurostd::xmatrix<int>() const;
  template JSON::object::operator aurostd::xmatrix<uint>() const;
  template JSON::object::operator aurostd::xmatrix<bool>() const;
  template JSON::object::operator aurostd::xmatrix<xcomplex<float>>() const;
  template JSON::object::operator aurostd::xmatrix<xcomplex<double>>() const;


  ///@brief allow to append to JSON::object_tpe::LIST
  void JSON::object::push_back(const JSON::object content){
    if (type == JSON::object_types::LIST) {
      std::shared_ptr <JSON::List> list_obj = std::static_pointer_cast<JSON::List>(obj);
      list_obj->push_back(content);
    }
    else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "push_back is just allowed for a JSON LIST", _VALUE_ILLEGAL_);
    }
  }

  ///@brief gives the size of JSON::object_tpe::LIST, JSON::object_tpe::DICTIONARY or JSON::object_tpe::STRING
  size_t JSON::object::size() {
    switch (type) {
      case object_types::DICTIONARY: {
        std::shared_ptr <JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
        return content->size();
      }
      case object_types::LIST: {
        std::shared_ptr <JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
        return content->size();
      }
      case object_types::STRING: {
        std::shared_ptr <std::string> content = std::static_pointer_cast<std::string>(obj);
        return content->size();
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "size is just allowed for a JSON DICTIONARY, LIST, or STRING", _VALUE_ILLEGAL_);
      }
    }
  }

  ///@brief checks if JSON::object_tpe::LIST, JSON::object_tpe::DICTIONARY or JSON::object_tpe::STRING is empty
  bool JSON::object::empty(){
    if ((type == object_types::DICTIONARY) || (type == object_types::LIST) || (type == object_types::STRING)) {
      return not bool(this);
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "empty is just allowed for a JSON DICTIONARY, LIST, or STRING", _VALUE_ILLEGAL_);
    }
  }
}
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
  static inline size_t range_find(const char *content_ptr, std::pair <size_t, size_t> border, char to_find) {
    return (char *) memchr(content_ptr + border.first, to_find, border.second - border.first + 1) - content_ptr;
  }

  /// @brief parse JSON string
  /// @param raw_content full json string
  /// @param border string borders
  /// @return c++ string
  /// @authors
  /// @mod{HE,20220924,created}
  std::string JSON::parse_string(const std::string &raw_content, std::pair <size_t, size_t> border)  {
    if (border.second == 0) {
      border.second = raw_content.size() - 1;
    }
    size_t last_escape_pos = border.first;
    size_t current_escape_pos = range_find(raw_content.c_str(), border, '\\');

    if (current_escape_pos > border.second) return raw_content.substr(border.first, border.second - border.first + 1);

    std::string result = "";
    while (current_escape_pos <= border.second) {
      result += raw_content.substr(last_escape_pos, current_escape_pos - last_escape_pos);
      switch (raw_content[current_escape_pos + 1]) {
        case '"': {
          result += '"';
        }
          break;
        case '\\': {
          result += '\\';
        }
          break;
        case '/': {
          result += '/';
        }
          break;
        case 'b': {
          result += '\b';
        }
          break;
        case 'f': {
          result += '\f';
        }
          break;
        case 'n': {
          result += '\n';
        }
          break;
        case 'r': {
          result += '\r';
        }
          break;
        case 't': {
          result += '\t';
        }
          break;
        case 'u': {
          if (current_escape_pos + 6 > raw_content.size())
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
          result += unescape_unicode(raw_content, current_escape_pos);
          current_escape_pos += 2;
        }
          break;
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
  std::string JSON::char_escape(const char16_t c){
    switch(c){
      default: return "";
      case('"'): return "\\\"";
      case('\\'): return "\\\\";
      case('/'): return "\\/";
      case('\b'): return "\\b";
      case('\f'): return "\\f";
      case('\n'): return "\\n";
      case('\r'): return "\\r";
      case('\t'): return "\\t";
    }
  }

  /// @brief prepare string for JSON output with or without unicode escapes
  /// @param raw string to escape
  /// @param unicode if true escape unicode charaters (default true)
  /// @return JSON string
  /// @authors
  /// @mod{HE,20221109,created}
  std::string JSON::escape(const std::string & raw, const bool unicode) {
    std::stringstream out;
    if (unicode) {
      std::u16string utf16 = std::wstring_convert < std::codecvt_utf8_utf16 < char16_t >, char16_t > {}.from_bytes(raw.data());
      for (char16_t c: utf16) {
        if (c < 0x80) {
          if (char_escape(c).empty()) out << (char) c;
          else out << char_escape(c);
          }
        else out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << c;
      }
    }
    else {
      for (char16_t c: raw) {
        if (char_escape(c).empty()) out << (char) c;
        else out << char_escape(c);
      }
    }
    return out.str();
  }

  /// @brief convert a unicode codepoint to a series of utf8 chars
  /// @param cp 32bit codepoint
  /// @return utf8 representation
  /// @authors
  /// @mod{HE,20221109,created}
  std::string JSON::char32_to_string(const char32_t cp)  {
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
  std::string JSON::unescape_unicode(const std::string & raw, size_t & pos) {
    pos += 2;
    char32_t cp = (char32_t) aurostd::string2utype<uint>(raw.substr(pos, 4), 16);

    if (cp < 0xd800 || cp > 0xdfff) return char32_to_string(cp);
    else if (cp > 0xdbff) throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
    else {
      if (raw[pos+4] == '\\' and raw[pos+5] == 'u'){
        pos += 6;
        char32_t trailing_cp = (char32_t) aurostd::string2utype<uint>(raw.substr(pos, 4), 16);
        if (trailing_cp < 0xdc00 || trailing_cp > 0xdfff) throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
        char32_t combo_cp = ((cp - 0xd800) << 10) + (trailing_cp - 0xdc00) + 0x10000;
        return char32_to_string(combo_cp);
      }
      else throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: undefined unicode character", _FILE_WRONG_FORMAT_);
    }
  }

  /// @brief strip whitespaces
  /// @param raw_content full json string
  /// @param border working border
  /// @return stripped raw string borders
  /// @authors
  /// @mod{HE,20220924,created}
  std::pair <size_t, size_t> JSON::find_strip(const std::string &raw_content, std::pair <size_t, size_t> border) {
    if (border.second == 0) border.second = raw_content.size() - 1;
    size_t start = raw_content.find_first_not_of(" \n\t\r\v\f", border.first);
    size_t end = raw_content.find_last_not_of(" \n\t\r\v\f", border.second);
    return {start, min(end, border.second)};
  }

  /// @brief find the border of a JSON string
  /// @param raw_content full json string
  /// @param border working border
  /// @return string borders
  /// @authors
  /// @mod{HE,20220924,created}
  std::pair <size_t, size_t> JSON::find_string(const std::string &raw_content, std::pair <size_t, size_t> border) {
    if (border.second == 0) border.second = raw_content.size() - 1;
    size_t start = range_find(raw_content.c_str(), border, '"');
    if (start >= border.second) return {std::string::npos, std::string::npos};
    size_t end = start;
    uint escape_check = 0;
    size_t escape_check_pos = 0;
    if (start == std::string::npos) return {start, end};
    do {
      end = range_find(raw_content.c_str(), {end + 1, border.second}, '"');
      escape_check_pos=end-1;
      escape_check = 0;
      while (raw_content[escape_check_pos] == '\\' and escape_check_pos >= start){
        escape_check++;
        escape_check_pos--;
      }
    } while ((escape_check%2!=0) && (end <= border.second));
    if (end > border.second)
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: string not closed", _FILE_WRONG_FORMAT_);
    return {start, end};
  }

  /// @brief find the border of a encapsulated by brackets
  /// @param raw_content full json string
  /// @param border working border
  /// @return bracket borders
  /// @authors
  /// @mod{HE,20220924,created}
  std::pair <size_t, size_t> JSON::find_bracket(const std::string &raw_content, char kind_open, std::pair <size_t, size_t> border) {
    char kind_close;
    switch (kind_open) {
      case '[':
        kind_close = ']';
        break;
      case '{':
        kind_close = '}';
        break;
      case '(':
        kind_close = ')';
        break;
      default:
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: bracket kind not allowed", _FILE_WRONG_FORMAT_);
    }
    size_t start = raw_content.find(kind_open, border.first);
    size_t end = start;
    if (start > border.second) return {std::string::npos, std::string::npos};
    size_t next_open=0;
    size_t next_close=0;
    std::pair <size_t, size_t> string_section;
    size_t open_count = 1;
    do {
//       cerr <<"end: " << end << " | "  << "open: " << open_count << " | current selection: " << raw_content.substr(start, end-start+1) << endl;
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
  JSON::object JSON::parse(const std::string &raw_content, std::pair <size_t, size_t> border) {
    if (border.second == 0) border.second = raw_content.size() - 1;
    object result;
    result.type = object_types::NONE;
    result.obj = nullptr;
    border = find_strip(raw_content, border);
    std::pair <size_t, size_t> section({0, 0});
    switch (raw_content[border.first]) {
      case '{': {
        if (raw_content[border.second] != '}')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: object not closed by '}'", _FILE_WRONG_FORMAT_);
        std::shared_ptr <JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
        border = {border.first + 1, border.second - 1};
        section = {0, 0};
        do {
          border = find_strip(raw_content, border);
//          cout  << border.first << " | "<< border.second << " - current: " << raw_content.substr(border.first, border.second-border.first+1) << endl;
          if (border.first >= border.second) break;
          section = find_string(raw_content, border);
          if (section.first == std::string::npos) break;
          std::string key_string = parse_string(raw_content, {section.first + 1, section.second - 1}); //cut " already
//          cout << "key: " << key_string << endl;
          border.first = raw_content.find(':', section.second) + 1;
          border = find_strip(raw_content, border);
//          cout  << border.first << " | "<< border.second << " - switch: " << raw_content.substr(border.first, border.second-border.first+1) << endl;

          switch (raw_content[border.first]) {
            case '"': {
              section = find_string(raw_content, border);
            }
              break;
            case '[': {
              section = find_bracket(raw_content, '[', border);
            }
              break;
            case '{': {
              section = find_bracket(raw_content, '{', border);
            }
              break;
            default: {
              section.first = border.first;
              section.second = min(border.second, range_find(raw_content.c_str(), border, ',') - 1);
            }
          }
          new_dictionary->insert({key_string, parse(raw_content, section)});
          border.first = range_find(raw_content.c_str(), {section.second, border.second}, ',');
          if (border.first == std::string::npos) border.first = border.second;
          else border.first++;
        } while (border.first < border.second);
        result.obj = new_dictionary;
        result.type = object_types::DICTIONARY;
      }
        break;
      case '[': {
        if (raw_content[border.second] != ']')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: list not closed by ']'", _FILE_WRONG_FORMAT_);
        std::shared_ptr <JSON::List> new_list = std::make_shared<JSON::List>();
        border = {border.first + 1, border.second - 1};
        section = {0, 0};
        do {
          border = find_strip(raw_content, border);
          switch (raw_content[border.first]) {
            case '"': {
              section = find_string(raw_content, border);
            }
              break;
            case '[': {
              section = find_bracket(raw_content, '[', border);
            }
              break;
            case '{': {
              section = find_bracket(raw_content, '{', border);
            }
              break;
            default: {
              section.first = border.first;
              section.second = min(border.second, raw_content.find(',', border.first) - 1);
            }
          }
          new_list->emplace_back(parse(raw_content, section));
          border.first = raw_content.find(',', section.second);
          if (border.first == std::string::npos) border.first = border.second;
          else border.first++;
        } while (section.second < border.second);
        result.obj = new_list;
        result.type = object_types::LIST;
      }
        break;
      case '"': {
        if (raw_content[border.second] != '"')
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON parsing failed: string not enclosed by '\\\"'\"", _FILE_WRONG_FORMAT_);
        std::shared_ptr <std::string> new_string = std::make_shared<std::string>();
        *new_string = parse_string(raw_content, {border.first + 1, border.second - 1});
        result.obj = new_string;
        result.type = object_types::STRING;
      }
        break;
      case 'n': {
        result.obj = nullptr;
        result.type = object_types::NONE;
      }
        break;
      case 't': {
        result.obj = nullptr;
        result.type = object_types::T;
      }
        break;
      case 'f': {
        result.obj = nullptr;
        result.type = object_types::F;
      }
        break;
      default: { //number
        if (raw_content.find('.', border.first) != std::string::npos || raw_content.find('E', border.first) != std::string::npos || raw_content.find('e', border.first) != std::string::npos) {
          std::shared_ptr<double> new_number = std::make_shared<double>();
          *new_number = std::strtod(raw_content.c_str() + border.first, nullptr);
          result.obj = new_number;
          result.type = object_types::FLOAT;
        } else {
          std::shared_ptr<long long> new_number = std::make_shared<long long>();
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
    std::string raw_content = aurostd::file2string(file_path);
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
  /// @param escape_unicode if true escape unicode characters
  /// @return converted string
  std::string JSON::toString(const object & root, const bool escape_unicode) {
    return root.toString(true, escape_unicode);
  }

  /// @brief save JSON::object to file
  /// @param root JSON::object to save
  /// @param file_path file path to save to
  /// @param escape_unicode if true escape unicode characters
  void JSON::saveFile(const object & root, const std::string & file_path, const bool escape_unicode) {
    aurostd::string2file(root.toString(true, escape_unicode), file_path);
  }

}

#endif // _AUROSTD_XPARSER_JSON_CPP_
