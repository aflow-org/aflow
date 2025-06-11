
#ifndef AUROSTD_XPARSER_JSON_TPP
#define AUROSTD_XPARSER_JSON_TPP

#include <algorithm>
#include <cstddef>
#include <deque>
#include <functional>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "aurostd.h"
#include "aurostd_xcomplex.h"
#include "aurostd_xerror.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xparser_json.h"
#include "aurostd_xvector.h"

namespace aurostd {

  /// @brief change this JSON::object to a JSON::object_types::INTEGER
  /// @authors
  /// @mod{HE,20220926,created}
  template <class utype> void JSON::object::fromNumber(const utype content) {
    if (std::is_integral<utype>::value) {
      this->type = object_types::INTEGER;
      const std::shared_ptr<long long> new_content = std::make_shared<long long>();
      *new_content = (long long) content;
      this->obj = new_content;
      return;
    } else if (std::is_floating_point<utype>::value) {
      if (std::isnan(content)) {
        this->type = object_types::NONE;
        this->obj = nullptr;
        return;
      }
      this->type = object_types::FLOAT;
      const std::shared_ptr<double> new_content = std::make_shared<double>();
      *new_content = static_cast<double>(content);
      this->obj = new_content;
      return;
    }
  }

  /// @brief change this JSON::object to a complex representation
  /// @authors
  /// @mod{HE,20221118,created}
  template <class utype> void JSON::object::fromComplex(const xcomplex<utype> &content) {
    const std::shared_ptr<JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
    this->obj = new_dictionary;
    this->type = object_types::DICTIONARY;
    new_dictionary->insert({"real", content.re});
    new_dictionary->insert({"imag", content.im});
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a xvector
  /// @authors
  /// @mod{HE,20220926,created}
  template <class utype> void JSON::object::fromXvector(const xvector<utype> &content) {
    const std::shared_ptr<JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (int idx = content.lrows; idx <= content.urows; idx++) {
      new_list->push_back(content[idx]);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a xmatrix
  /// @authors
  /// @mod{HE,20220926,created}
  template <class utype> void JSON::object::fromXmatrix(const xmatrix<utype> &content) {
    const std::shared_ptr<JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (int idx = content.lrows; idx <= content.urows; idx++) {
      new_list->push_back(content(idx));
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a vector
  /// @authors
  /// @mod{HE,20220926,created}
  template <class utype> void JSON::object::fromVector(const std::vector<utype> &content) {
    const std::shared_ptr<JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (utype entry : content) {
      new_list->push_back(entry);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on a deque
  /// @authors
  /// @mod{HE,20240222,created}
  template <class utype> void JSON::object::fromDeque(const std::deque<utype> &content) {
    const std::shared_ptr<JSON::List> new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = object_types::LIST;
    for (utype entry : content) {
      new_list->push_back(entry);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::DICTIONARY based on a map
  /// @authors
  /// @mod{HE,20220926,created}
  template <class utype> void JSON::object::fromMap(const std::map<std::string, utype> &content) {
    const std::shared_ptr<JSON::Dictionary> new_dictionary = std::make_shared<JSON::Dictionary>();
    this->obj = new_dictionary;
    this->type = object_types::DICTIONARY;
    for (const std::pair<std::string, utype> entry : content) {
      new_dictionary->insert({entry.first, entry.second});
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::LIST based on list like types
  /// See @c MockSerializableListDict for tested types.
  /// @authors
  /// @mod{ST,20250315,created}
  template <class T, JSON::enable_list_like<T>> void JSON::object::fromListLike(const T &content) {
    const auto new_list = std::make_shared<JSON::List>();
    this->obj = new_list;
    this->type = JSON::object_types::LIST;
    for (const auto &entry : content) {
      new_list->push_back(entry);
    }
  }

  /// @brief change this JSON::object to a JSON::object_types::DICTIONARY based on a map
  /// SFINAE is volatile here, the JSON::enable_dict_like type trait does not fully protect the duck typing
  /// This supports types that are strictly map-like with string keys, such as map and unordered_map, but support
  /// for others may not be guaranteed. Test new types. See @c MockSerializableListDict for tested types.
  /// @authors
  /// @mod{ST,20250315,created}
  template <class T, JSON::enable_dict_like<T>> void JSON::object::fromDictLike(const T &content) {
    const auto new_dictionary = std::make_shared<JSON::Dictionary>();
    this->obj = new_dictionary;
    this->type = JSON::object_types::DICTIONARY;
    for (const auto &[key, value] : content) {
      new_dictionary->insert({key, value});
    }
  }

  /// @brief converting constructor: set the content of this JSON::object based on a number
  template <typename utype, JSON::enable_arithmetic<utype>> JSON::object::object(const utype content) {
    this->fromNumber(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a number
  template <class utype> JSON::object::object(const xcomplex<utype> &content) {
    this->fromComplex(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a xvector
  template <class utype> JSON::object::object(const xvector<utype> &content) {
    this->fromXvector(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a xmatrix
  template <class utype> JSON::object::object(const xmatrix<utype> &content) {
    this->fromXmatrix(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a vector
  template <class utype> JSON::object::object(const vector<utype> &content) {
    this->fromVector(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a deque
  template <class utype> JSON::object::object(const deque<utype> &content) {
    this->fromDeque(content);
  }

  /// @brief converting constructor: set the content of this JSON::object based on a map
  template <class utype> JSON::object::object(const std::map<std::string, utype> &content) {
    this->fromMap(content);
  }

  template <class T, JSON::enable_list_like<T>> JSON::object::object(const T &content) {
    this->fromListLike(content);
  }

  template <class T, JSON::enable_dict_like<T>> JSON::object::object(const T &content) {
    this->fromDictLike(content);
  }

  ///@brief assignment operator for numbers
  template <class utype, JSON::enable_arithmetic<utype>> JSON::object &JSON::object::operator=(const utype content) {
    this->fromNumber(content);
    return *this;
  }

  ///@brief assignment operator for xcomplex
  template <class utype> JSON::object &JSON::object::operator=(const xcomplex<utype> &content) {
    fromComplex(content);
    return *this;
  }

  ///@brief assignment operator for xvector
  template <class utype> JSON::object &JSON::object::operator=(const xvector<utype> &content) {
    fromXvector(content);
    return *this;
  }

  ///@brief assignment operator for xmatrix
  template <class utype> JSON::object &JSON::object::operator=(const xmatrix<utype> &content) {
    fromXmatrix(content);
    return *this;
  }

  ///@brief assignment operator for vector
  template <class utype> JSON::object &JSON::object::operator=(const std::vector<utype> &content) {
    fromVector(content);
    return *this;
  }

  ///@brief assignment operator for deque
  template <class utype> JSON::object &JSON::object::operator=(const std::deque<utype> &content) {
    fromDeque(content);
    return *this;
  }

  ///@brief assignment operator for map
  template <class utype> JSON::object &JSON::object::operator=(const std::map<std::string, utype> &content) {
    fromMap(content);
    return *this;
  }

  template <class utype> JSON::object &JSON::object::operator=(const JsonSerializable<utype> &content) {
    return operator=(content.serialize());
  }

  template <class T, JSON::enable_list_like<T>> JSON::object &JSON::object::operator=(const T &content) {
    this->fromListLike(content);
    return *this;
  }

  template <class T, JSON::enable_dict_like<T>> JSON::object &JSON::object::operator=(const T &content) {
    this->fromDictLike(content);
    return *this;
  }

  ///@brief conversion function for int
  template <class utype> JSON::object::operator xcomplex<utype>() const {
    if (type != object_types::DICTIONARY) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xcomplex conversion failed: element is not a DICTIONARY: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    return xcomplex<utype>(static_cast<utype>(content->operator[]("real")), static_cast<utype>(content->operator[]("imag")));
  }

  ///@brief conversion function for vector
  template <class utype> JSON::object::operator std::vector<utype>() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    std::vector<utype> result;
    const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) {
      return result;
    }
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (const JSON::object &so : *content) {
        result.emplace_back(static_cast<utype>(so));
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (const JSON::object &so : *content) {
        if (so.type == object_types::DICTIONARY) {
          result.push_back(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING) {
          result.push_back(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) {
          result.push_back(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else {
      // fallback to generic
      return toListLike<std::vector<utype>>();
    }
    return result;
  }

  ///@brief conversion function for deque
  template <class utype> JSON::object::operator std::deque<utype>() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON deque conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    std::deque<utype> result;
    const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) {
      return result;
    }
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (const JSON::object &so : *content) {
        result.emplace_back(static_cast<utype>(so));
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (const JSON::object &so : *content) {
        if (so.type == object_types::DICTIONARY) {
          result.push_back(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING) {
          result.push_back((utype) so);
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON deque conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) {
          result.push_back(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON deque conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else {
      // fallback to generic
      return toListLike<std::deque<utype>>();
    }
    return result;
  }

  ///@brief conversion function for set with type conversion
  template <class utype> JSON::object::operator std::set<utype>() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON set conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    std::set<utype> result;
    const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) {
      return result;
    }
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (const JSON::object &so : *content) {
        result.emplace(static_cast<utype>(so));
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (const JSON::object &so : *content) {
        if (so.type == object_types::DICTIONARY) {
          result.emplace(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING) {
          result.emplace(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) {
          result.emplace(static_cast<utype>(so));
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON vector conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
      }
    } else {
      // fallback to generic
      for (const auto &entry : *content) {
        result.emplace(entry);
      }
    }
    return result;
  }

  ///@brief conversion function for map with type conversion
  template <class utype> JSON::object::operator std::map<std::string, utype>() const {
    if (type != object_types::DICTIONARY) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a DICTIONARY: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    std::map<std::string, utype> result;
    const std::shared_ptr<JSON::Dictionary> content = std::static_pointer_cast<JSON::Dictionary>(obj);
    if (content->empty()) {
      return result;
    }
    if (typeid(utype) == typeid(std::string) || typeid(utype) == typeid(bool)) {
      for (auto &[key, value] : *content) {
        result.insert({key, static_cast<utype>(value)});
      }
    } else if (std::is_floating_point<utype>::value) {
      for (auto &[key, value] : *content) {
        if (value.type > object_types::STRING) {
          result.insert({key, static_cast<utype>(value)});
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + static_cast<std::string>(value), _VALUE_ILLEGAL_);
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (auto &[key, value] : *content) {
        if (value.type > object_types::STRING && value.type < object_types::NONE) {
          result.insert({key, static_cast<utype>(value)});
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON std::map conversion failed: element is not a number: " + static_cast<std::string>(value), _VALUE_ILLEGAL_);
        }
      }
    } else {
      // fallback to generic
      return toDictLike<std::map<std::string, utype>>();
    }
    return result;
  }

    ///@brief conversion function for xvector
  template <class utype> JSON::object::operator aurostd::xvector<utype>() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xvector conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    const aurostd::xvector<utype> result(content->size(), 1);
    if (content->empty()) {
      return result;
    }
    size_t idx = 1;
    if (typeid(utype) == typeid(bool)) {
      for (const JSON::object &so : *content) {
        result[idx] = static_cast<utype>(so);
        idx++;
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (const JSON::object &so : *content) {
        if (so.type == object_types::DICTIONARY) {
          result[idx] = static_cast<utype>(so);
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
        idx++;
      }
    } else if (std::is_floating_point<utype>::value) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING) {
          result[idx] = static_cast<utype>(so);
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
        idx++;
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (const JSON::object &so : *content) {
        if (so.type > object_types::STRING && so.type < object_types::NONE) {
          result[idx] = static_cast<utype>(so);
        } else {
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a number: " + static_cast<std::string>(so), _VALUE_ILLEGAL_);
        }
        idx++;
      }
    }
    return result;
  }

  ///@brief conversion function for xmatrix
  template <class utype> JSON::object::operator aurostd::xmatrix<utype>() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) {
      return aurostd::xmatrix<utype>();
    }

    // size scan
    size_t rows = 0;
    size_t cols = 0;
    std::vector<size_t> cols_v;
    for (const JSON::object &so : *content) {
      if (so.type != object_types::LIST) {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
      }
      const std::shared_ptr<JSON::List> row = std::static_pointer_cast<JSON::List>(so.obj);
      rows++;
      cols_v.emplace_back(row->size());
    }
    if (std::adjacent_find(cols_v.begin(), cols_v.end(), std::not_equal_to<size_t>()) == cols_v.end()) {
      cols = cols_v[0];
    } else {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: sub elements have different length", _VALUE_ILLEGAL_);
    }

    const aurostd::xmatrix<utype> result(rows, cols);
    if (typeid(utype) == typeid(bool)) {
      for (int r = result.lrows; r <= result.urows; r++) {
        const std::shared_ptr<JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          result[r][c] = static_cast<utype>(row->operator[](c - 1));
        }
      }
    } else if ((typeid(utype) == typeid(aurostd::xcomplex<double>)) || (typeid(utype) == typeid(aurostd::xcomplex<float>))) {
      for (int r = result.lrows; r <= result.urows; r++) {
        const std::shared_ptr<JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type == object_types::DICTIONARY) {
            result[r][c] = (utype) row->operator[](c - 1);
          } else {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a complex number: " + static_cast<std::string>(row->operator[](c - 1)), _VALUE_ILLEGAL_);
          }
        }
      }
    } else if (std::is_floating_point<utype>::value) {
      for (int r = result.lrows; r <= result.urows; r++) {
        const std::shared_ptr<JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING) {
            result[r][c] = static_cast<utype>(row->operator[](c - 1));
          } else {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + static_cast<string>(row->operator[](c - 1)), _VALUE_ILLEGAL_);
          }
        }
      }
    } else if (std::numeric_limits<utype>::is_integer) {
      for (int r = result.lrows; r <= result.urows; r++) {
        const std::shared_ptr<JSON::List> row = std::static_pointer_cast<JSON::List>(content->operator[](r - 1).obj);
        for (int c = result.lcols; c <= result.ucols; c++) {
          if (row->operator[](c - 1).type > object_types::STRING && row->operator[](c - 1).type < object_types::NONE) {
            result[r][c] = static_cast<utype>(row->operator[](c - 1));
          } else {
            throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON xmatrix conversion failed: elements is not a number: " + static_cast<string>(row->operator[](c - 1)), _VALUE_ILLEGAL_);
          }
        }
      }
    }
    return result;
  }

  /// @brief conversion function for list-like types (e.g. vector, deque, list)
  /// See @c MockSerializableListDict for tested types.
  /// @authors
  /// @mod{ST,20250315,created}
  template <class T, JSON::enable_list_like<T>> T JSON::object::toListLike() const {
    if (type != object_types::LIST) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a LIST: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    T result;
    const auto content = std::static_pointer_cast<JSON::List>(obj);
    if (content->empty()) {
      return result;
    }
    for (const auto &entry : *content) {
      result.emplace_back(entry);
    }
    return result;
  }

  /// @brief conversion function for dict-like types (e.g. map, hashmap) where the key is a string
  /// See @c MockSerializableListDict for tested types.
  /// @authors
  /// @mod{ST,20250315,created}
  template <class T, JSON::enable_dict_like<T>> T JSON::object::toDictLike() const {
    if (type != object_types::DICTIONARY) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "JSON conversion failed: element is not a DICTIONARY: " + static_cast<std::string>(*this), _VALUE_ILLEGAL_);
    }
    T result;
    const auto content = std::static_pointer_cast<JSON::Dictionary>(obj);
    if (content->empty()) {
      return result;
    }
    for (const auto &[key, value] : *content) {
      result.emplace(key, value);
    }
    return result;
  }

  template <class T, JSON::enable_list_like<T>> JSON::object::operator T() const {
    return toListLike<T>();
  }

  template <class T, JSON::enable_dict_like<T>> JSON::object::operator T() const {
    return toDictLike<T>();
  }

  template <class T, JSON::enable_serializable<T>> JSON::object::operator T() const {
    T t;
    JsonSerializable<T> &tr = t;
    // dodge protected access because JsonSerializable is friend
    return tr.deserialize(*this);
  }

  template <class utype> typename std::list<utype>::iterator JSON::object::find_list(const utype &key) const {
    switch (type) {
      case object_types::LIST: {
        const std::shared_ptr<JSON::List> content = std::static_pointer_cast<JSON::List>(obj);
        return std::find(content->begin(), content->end(), key);
      }
      default: {
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "find_list is just allowed for a JSON LIST", _VALUE_ILLEGAL_);
      }
    }
  }

} // namespace aurostd

#endif // AUROSTD_XPARSER_JSON_TPP
