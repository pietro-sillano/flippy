#ifndef FLIPPY_CUSTOM_CONCEPTS_HPP
#define FLIPPY_CUSTOM_CONCEPTS_HPP
#include <concepts>

namespace fp{
/**
 * Here we implement the concepts of a floating.
 * This concept is equivalent to [std::floating_point](https://en.cppreference.com/w/cpp/concepts/floating_point).
 * However, not all `c++20` compilers implement these features yet.
 * In particular some apple clang versions that are technically `c++20` capable.
 *
 * @note For more informaton on what `c++` concepts are in general see [this](https://en.cppreference.com/w/cpp/language/constraints) cpp reference page.
 */
template<class T> concept floating_point_number = std::is_floating_point_v<T>;

/**
 * Here we implement the concepts of an integer number.
 * This concept is equivalent to [std::integral](https://en.cppreference.com/w/cpp/concepts/integral).
 * However, not all `c++20` compilers implement these features yet.
 * In particular some apple clang versions that are technically `c++20` capable.
 *
 * @note For more informaton on what `c++` concepts are in general see [this](https://en.cppreference.com/w/cpp/language/constraints) cpp reference page.
 */
template<class T> concept integer_number = std::is_integral_v<T>;
}


#endif //FLIPPY_CUSTOM_CONCEPTS_HPP
