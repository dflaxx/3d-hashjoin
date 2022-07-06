#pragma once

#include <string_view>	

// Macros for printing variables, variable names, types (poor man's reflection)
// Source: Magnus

#define COUT_LOG std::cout << __FILE__ << ":" << __LINE__ << ": "

#define STRINGIFY(X) #X
#define EXPAND_AND_STRINGIFY(X) STRINGIFY(X)
#define PRINT_REFLECTION(X) (#X) << ": " << (X)
#define TRACE_AND_EXECUTE(X) COUT_LOG << EXPAND_AND_STRINGIFY(X) << std::endl; X;

// Get typename as string in human-readable form
// usage:  int x;
//         std::cout << type_name<decltype(x)>() << std::endl;  // prints 'int'
// source: https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c/
template <typename T>
constexpr auto type_name() noexcept {
  std::string_view name = "Error: unsupported compiler", prefix, suffix;
#ifdef __clang__
  name = __PRETTY_FUNCTION__;
  prefix = "auto type_name() [T = ";
  suffix = "]";
#elif defined(__GNUC__)
  name = __PRETTY_FUNCTION__;
  prefix = "constexpr auto type_name() [with T = ";
  suffix = "]";
#elif defined(_MSC_VER)
  name = __FUNCSIG__;
  prefix = "auto __cdecl type_name<";
  suffix = ">(void) noexcept";
#endif
  name.remove_prefix(prefix.size());
  name.remove_suffix(suffix.size());
  return name;
}

