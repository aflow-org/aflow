
#ifndef AUROSTD_TYPE_TRAITS_TPP
#define AUROSTD_TYPE_TRAITS_TPP

#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace aurostd {

// Metaprogramming helpers

// pair

  template <typename T> struct is_pair_like : std::false_type {};

  template <typename T, typename U> struct is_pair_like<std::pair<T, U>> : std::true_type {};

  template <typename T> inline constexpr bool is_pair_like_v = is_pair_like<T>::value;

  template <typename T> using pair_first_t = std::tuple_element_t<0, T>;

  template <typename T> using pair_second_t = std::tuple_element_t<1, T>;

/// Check if the first type in pair T is the given type F1
  template <typename T, typename F1> struct is_pair_with_first : std::conjunction<is_pair_like<T>, std::is_same<std::decay_t<pair_first_t<T>>, std::decay_t<F1>>> {};

  template <typename T, typename F1> inline constexpr bool is_pair_with_first_v = is_pair_with_first<T, F1>::value;

/// Check if the first type in pair T is std::string
  template <typename T> struct is_string_first_pair : is_pair_with_first<T, std::string> {};

  template <typename T> inline constexpr bool is_string_first_pair_v = is_string_first_pair<T>::value;

// tuple

  template <typename T> struct is_tuple_like : std::false_type {};

  template <typename... T> struct is_tuple_like<std::tuple<T...>> : std::true_type {};

  template <typename T> inline constexpr bool is_tuple_like_v = is_tuple_like<T>::value;

  template <size_t N, typename F, typename T> struct tuple_element_is : std::conjunction<is_tuple_like<T>, std::is_same<std::tuple_element_t<N, std::decay_t<T>>, F>> {};

  template <size_t N, typename F, typename T> inline constexpr bool tuple_element_is_v = tuple_element_is<N, F, T>::value;

// iterator

/// Check for type that has an iterator
  template <typename T, typename = void> struct has_iterator : std::false_type {};

/// Check for type that has an iterator
  template <typename T> struct has_iterator<T, std::void_t<typename std::decay_t<T>::iterator>> : std::true_type {};

  template <typename T> inline constexpr bool has_iterator_v = has_iterator<T>::value;

// list-like iterables

/// Check for type that has an iterator over a type that is not a pair
  template <typename T, typename = void> struct is_list_like : std::false_type {};

/// Check for type that has an iterator over a type that is not a pair
  template <typename T> struct is_list_like<T, std::void_t<typename std::decay_t<T>::iterator>> : std::negation<is_pair_like<typename std::decay_t<T>::value_type>> {};

  template <typename T> inline constexpr bool is_list_like_v = is_list_like<T>::value;

// map-like iterables

/// Check for type that has an iterator over a type that is a pair
  template <typename T, typename = void> struct is_map_like : std::false_type {};

/// Check for type that has an iterator over a type that is a pair
  template <typename T> struct is_map_like<T, std::void_t<typename std::decay_t<T>::iterator>> : is_pair_like<typename std::decay_t<T>::value_type> {};

  template <typename T> inline constexpr bool is_map_like_v = is_map_like<T>::value;

// dict-like iterables

/// Check for type that has an iterator over a type that is a pair and whose first type is a string
  template <typename T, typename = void> struct is_dict_like : std::false_type {};

/// Check for type that has an iterator over a type that is a pair and whose first type is a string
  template <typename T>
  struct is_dict_like<T, std::void_t<typename std::decay_t<T>::iterator>> : std::conjunction<is_pair_like<typename std::decay_t<T>::value_type>, is_string_first_pair<typename std::decay_t<T>::value_type>> {};

  template <typename T> inline constexpr bool is_dict_like_v = is_dict_like<T>::value;

} // namespace aurostd

#endif // AUROSTD_TYPE_TRAITS_TPP
