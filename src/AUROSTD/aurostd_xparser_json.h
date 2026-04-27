
#ifndef AUROSTD_XPARSER_JSON_H
#define AUROSTD_XPARSER_JSON_H

#include <cmath>
#include <cstddef>
#include <deque>
#include <iosfwd>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "aurostd_type_traits.tpp"
#include "aurostd_xcomplex.h"
#include "aurostd_xfile.h"
#include "aurostd_xmatrix.h"
#include "aurostd_xvector.h"

template <typename T> class JsonSerializable;

namespace aurostd {
  // JSON namespace for reading and writing //HE20221109
  namespace JSON {

    template <typename T> using is_string_like = std::is_same<std::decay_t<T>, std::string>;
    template <typename T> using is_json_list_like = std::conjunction<aurostd::is_list_like<T>, std::negation<is_string_like<T>>>;
    template <typename T> using enable_list_like = std::enable_if_t<is_json_list_like<T>::value, bool>;
    template <typename T> using enable_dict_like = std::enable_if_t<aurostd::is_dict_like_v<T>, bool>;
    template <typename T> using enable_serializable = std::enable_if_t<std::is_base_of_v<JsonSerializable<T>, T>, bool>;
    template <typename T> using enable_arithmetic = std::enable_if_t<std::is_arithmetic_v<T>, bool>;

    typedef std::vector<object> List;                 ///< shortcut for JSON::object_types::LIST
    typedef std::map<std::string, object> Dictionary; ///< shortcut for JSON::object_types::DICTIONARY

    enum class object_types { // order is important!
      DICTIONARY,
      LIST,
      STRING,
      FLOAT,
      INTEGER,
      T,
      F,
      NONE
    };

    struct object {
      object_types type = object_types::NONE;
      std::shared_ptr<void> obj = nullptr;

      // operators
      JSON::object& operator[](size_t index);
      const JSON::object& operator[](size_t index) const;
      JSON::object& operator[](const std::string& key);
      const JSON::object& operator[](const std::string& key) const;
      JSON::object& operator[](const char* key);
      const JSON::object& operator[](const char* key) const;
      JSON::object& operator=(const char* content); // for literal strings
      JSON::object& operator=(const std::string& content);
      JSON::object& operator=(bool content);
      JSON::object& operator=(std::nullptr_t content);
      JSON::object& operator=(const List& content);
      JSON::object& operator=(const Dictionary& content);
      template <class utype, enable_arithmetic<utype> = true> JSON::object& operator=(utype content);
      template <class utype> JSON::object& operator=(const xcomplex<utype>& content);
      template <class utype> JSON::object& operator=(const std::vector<utype>& content);
      template <class utype> JSON::object& operator=(const std::deque<utype>& content);
      template <class utype> JSON::object& operator=(const std::map<std::string, utype>& content);
      template <class utype> JSON::object& operator=(const xvector<utype>& content);
      template <class utype> JSON::object& operator=(const xmatrix<utype>& content);
      template <class utype> JSON::object& operator=(const JsonSerializable<utype>& content);
      template <class T, enable_list_like<T> = true> JSON::object& operator=(const T& content);
      template <class T, enable_dict_like<T> = true> JSON::object& operator=(const T& content);

      // converting constructors
      object() = default;
      object(const char* content);
      object(const std::string& content);
      object(bool content);
      object(std::nullptr_t content);
      object(object_types create_type);
      object(const List& content);
      object(const Dictionary& content);
      template <typename utype, enable_arithmetic<utype> = true> object(utype content);
      template <typename utype> object(const xcomplex<utype>& content);
      template <typename utype> object(const std::vector<utype>& content);
      template <typename utype> object(const std::deque<utype>& content);
      template <typename utype> object(const std::map<std::string, utype>& content);
      template <typename utype> object(const xvector<utype>& content);
      template <typename utype> object(const xmatrix<utype>& content);
      template <typename T> object(const JsonSerializable<T>& content) : JSON::object(std::move(content.serialize())) {}
      template <typename T, enable_list_like<T> = true> object(const T& content);
      template <typename T, enable_dict_like<T> = true> object(const T& content);

      // conversion functions
      explicit operator bool() const;
      explicit operator std::string() const;
      explicit operator double() const;
      explicit operator long double() const;
      explicit operator float() const;
      explicit operator long long() const;
      explicit operator long() const;
      explicit operator int() const;
      explicit operator char() const;
      explicit operator unsigned long long() const;
      explicit operator unsigned long() const;
      explicit operator unsigned int() const;
      operator List() const;
      operator Dictionary() const;
      template <class utype> operator xcomplex<utype>() const;
      template <class utype> operator std::vector<utype>() const;
      template <class utype> operator std::deque<utype>() const;
      template <class utype> operator std::set<utype>() const;
      template <class utype> operator std::map<std::string, utype>() const;
      template <class utype> operator aurostd::xvector<utype>() const;
      template <class utype> operator aurostd::xmatrix<utype>() const;
      template <class T, enable_list_like<T> = true> explicit operator T() const;
      template <class T, enable_dict_like<T> = true> explicit operator T() const;
      template <class T, enable_serializable<T> = true> explicit operator T() const;

      // type specific functions
      void push_back(const JSON::object& content) const;
      void join(const JSON::object& content) const;
      [[nodiscard]] size_t size() const;
      [[nodiscard]] bool empty() const;
      [[nodiscard]] size_t count(const std::string&) const;
      [[nodiscard]] std::map<std::string, JSON::object>::iterator find(const std::string& key) const;
      template <class utype> typename std::list<utype>::iterator find_list(const utype& key) const;
      [[nodiscard]] std::map<std::string, JSON::object>::iterator end() const;
      [[nodiscard]] std::map<std::string, JSON::object>::iterator begin() const;
      [[nodiscard]] std::vector<std::string> keys() const;

      // conversion helper
      void fromString(const std::string& content);
      void fromList(const List& content);
      void fromDictionary(const Dictionary& content);
      template <class utype> void fromNumber(utype content);
      template <class utype> void fromComplex(const xcomplex<utype>& content);
      template <class utype> void fromVector(const std::vector<utype>& content);
      template <class utype> void fromDeque(const std::deque<utype>& content);
      template <class utype> void fromMap(const std::map<std::string, utype>& content);
      template <class utype> void fromXvector(const xvector<utype>& content);
      template <class utype> void fromXmatrix(const xmatrix<utype>& content);
      template <class T, enable_list_like<T> = true> void fromListLike(const T& content);
      template <class T, enable_dict_like<T> = true> void fromDictLike(const T& content);
      template <class T, enable_list_like<T> = true> T toListLike() const;
      template <class T, enable_dict_like<T> = true> T toDictLike() const;

      [[nodiscard]] std::string toString(bool json_format = true, bool escape_unicode = true) const;
      void saveFile(const std::string& file_path, compression_type ct = compression_type::None, bool escape_unicode = true) const;
    };

    // unicode helper function
    std::string unescape_unicode(const std::string& raw, size_t& pos);
    std::string escape(const std::string& raw, bool unicode = true);
    std::string char32_to_string(char32_t cp);
    std::string char_escape(char32_t c);

    // basic functions
    object loadFile(const std::string& file_path);
    object loadString(const std::string& content);
    void saveFile(const object& root, const std::string& file_path, compression_type ct = compression_type::None, bool escape_unicode = true);
    std::string toString(const object& root, bool escape_unicode = false);

    // navigation functions
    std::pair<size_t, size_t> find_string(const std::string& raw_content, std::pair<size_t, size_t> border = {0, 0});
    std::pair<size_t, size_t> find_bracket(const std::string& raw_content, char kind_open, std::pair<size_t, size_t> border = {0, 0});
    std::pair<size_t, size_t> find_strip(const std::string& raw_content, std::pair<size_t, size_t> border = {0, 0});

    // parser core
    object parse(const std::string& raw_content, std::pair<size_t, size_t> border = {0, 0});
    std::string parse_string(const std::string& raw_content, std::pair<size_t, size_t> border = {0, 0});

  }; // namespace JSON

  // insertion operator: enable easy interaction with cout
  std::ostream& operator<<(std::ostream& os, const JSON::object& jo);
} // namespace aurostd

#include "aurostd_xparser_json.tpp" // NOLINT / template implementation

#endif // AUROSTD_XPARSER_JSON_H
