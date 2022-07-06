/*
 * This implementation uses C++20 constraints and concepts
 * to enforce compile-time requirements on class templates and template arguments.
 *
 * Used for the hash table implementations in ht_chaining.hh and ht_nested.hh,
 * and the physical algebra operators in algebra.hh.
 */
#pragma once

#include "util/standard_includes.hh"

#if __cplusplus >= 202002L
#include <concepts>

// Printable: has operator<<(std::ostream&)
template <typename T>
concept Printable = requires(std::ostream& os, T a) {
  os << a;
};

// Hash function: has static eval(const input_t*) -> output_t
template <typename T>
concept alg_hashfun_c = requires(T t) {
  typename T::input_t;
  typename T::output_t;
  { T::eval(static_cast<const typename T::input_t*>(nullptr)) }
    -> std::same_as<typename T::output_t>;
};

// unary predicate (static): has static eval(const input_t*) -> bool
template <typename T>
concept alg_predicate_c = requires(T t) {
  typename T::input_t;
  { T::eval(static_cast<const typename T::input_t*>(nullptr)) }
    -> std::same_as<bool>;
};

// unary predicate (dynamic): has operator()(bool input_t*) -> bool
template <typename T>
concept alg_dyn_predicate_c = requires(T t) {
  typename T::input_t;
  // must have operator()
  { t(static_cast<const typename T::input_t*>(nullptr)) }
    -> std::same_as<bool>;
};

// binary predicate (static): has static eval(const input_t*) -> bool
// (used for key comparison function and join predicate (probe))
template <typename T>
concept alg_binary_predicate_c = requires(T t) {
  typename T::left_t;
  typename T::right_t;
  { T::eval(static_cast<const typename T::left_t*>(nullptr),
            static_cast<const typename T::right_t*>(nullptr)) }
    -> std::same_as<bool>;
};

// concatenation function to concatenate two database tuples:
// has function static eval(left_t*, right_t*) -> output_t
template <typename T>
concept alg_concatfun_c = requires(T t) {
  typename T::left_t;
  typename T::right_t;
  typename T::output_t;
  { T::eval(static_cast<typename T::left_t*>(nullptr),
            static_cast<typename T::right_t*>(nullptr)) }
    -> std::same_as<typename T::output_t>;
};

// unnest function for 3D hash join on nested hash table:
template <typename T>
concept alg_unnestfun_c = requires(T t) {
  typename T::input_t;  // a nested tuple type
  typename T::output_t; // an unnested tuple type
  typename T::MainNode; // from HtNested
  typename T::data_t;   // from HtNested
  { T::eval_left(static_cast<typename T::output_t*>(nullptr),
                 static_cast<typename T::input_t*>(nullptr)) }
    -> std::same_as<void>;
  { T::eval_right(static_cast<typename T::output_t*>(nullptr),
                  static_cast<typename T::input_t*>(nullptr),
                  static_cast<const typename T::data_t*>(nullptr)) }
    -> std::same_as<void>;
  { T::getMainNode(static_cast<typename T::input_t*>(nullptr)) }
    -> std::same_as<const typename T::MainNode*>;
};
#else
#error Compiling with C++ < 20, concepts not supported.
#endif  // __cplusplus >= 202002L
