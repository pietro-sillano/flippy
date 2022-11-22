#ifndef FLIPPY_CUSTOM_CONCEPTS_HPP
#define FLIPPY_CUSTOM_CONCEPTS_HPP
#include <concepts>

namespace fp{
/**
 * @defgroup myConcepts Special concepts
 * @brief definitions of concepts that restrict template parameters in  flippy's classes.
 * @note For more information on what `c++` concepts are in general see [this](https://en.cppreference.com/w/cpp/language/constraints) cpp reference page.
 *@{
 */

/**
 * @brief Here we implement the concepts of a floating point number.
 *
 * This concept is equivalent to [std::floating_point](https://en.cppreference.com/w/cpp/concepts/floating_point).
 * However, not all `c++20` compilers implement these features yet.
 * In particular some apple clang versions that are technically `c++20` capable.
 * @tparam T This concept requires the type T to be a floating point number.
 */
template<class T> concept floating_point_number = std::is_floating_point_v<T>;

/**
 * @brief Here we implement the concepts of a positive integer number that is used throughout the code for indexing.
 *
 * @tparam T This concept requires the type T to be unsigned and an integral type.
 */
template<class T> concept indexing_number = std::is_unsigned_v<T> && std::is_integral_v<T>;
/**@}*/
}


#endif //FLIPPY_CUSTOM_CONCEPTS_HPP
