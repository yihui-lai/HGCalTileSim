/**
 * @file    VectorUtils.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Template functions for common vector operations.
 */

#ifndef USERUTILS_COMMON_STLUTILS_VECTORUTILS_HPP
#define USERUTILS_COMMON_STLUTILS_VECTORUTILS_HPP
#include <vector>
#include <set>
#include <algorithm>

namespace usr {

/**
 * @defgroup VecUtils Vector enhancements
 * @brief    std::vector enhancement functions.
 * @ingroup  Common
 * @details  A collection of STL recipes for performing vector operations, many
 * of which consists of full range parsing of vectors. Another important factor
 * is the variadic template interface for creating vectors, something which
 * will enable the use of variadic template interface for other function.
 * @{
 */
template<typename T>
std::vector<T> RemoveDuplicate( const std::vector<T>& dupvec );

template<typename OBJTYPE, typename ARGTYPE, typename ... VARTYPE>
std::vector<OBJTYPE> MakeVector( const ARGTYPE& first, const VARTYPE& ... arg );

template<typename OBJTYPE, typename ARGTYPE>
std::vector<OBJTYPE> MakeVector( const ARGTYPE& first );

template<typename OBJTYPE>
void ClearValue( std::vector<OBJTYPE>& vec, const OBJTYPE& x );

template<typename OBJTYPE>
bool FindValue( std::vector<OBJTYPE>& vec, const OBJTYPE& x );
/** @} */

/*-----------------------------------------------------------------------------
 *  Template Implementations
   --------------------------------------------------------------------------*/
/**
 * @brief Removing duplicates from a vector while maintaining vector order
 * regardless of whether the vector is sorted or not.
 */
template<typename T>
std::vector<T>
RemoveDuplicate( const std::vector<T>& dupvec )
{
  std::set<T> unique_set;
  std::vector<T> ans;

  for( const auto& item : dupvec ){
    if( unique_set.insert( item ).second ){
      ans.push_back( item );
    }
  }

  return ans;
}

/**
 * @brief Constructing a vector of same-typed inputs of arbitrary number of
 *        function inputs
 *
 * Basically, this functions allow for the function call like:
 * `MakeVector<int>(12,3,4,5,6,7,8,9)`, which is turn could enable functions
 * to be declared as:
 * ```cpp
 * f( int reqired, TYPES... args){ return f(required, MakeVector(args...)); }
 * ```
 * allowing the user to input arbitrary many arguments for a given function,
 * particularly useful if your functions has many defaulted arguments that
 * the user might not want to keep track of ordering.
 */
template<typename OBJTYPE, typename ARGTYPE, typename ... VARTYPE>
std::vector<OBJTYPE>
MakeVector( const ARGTYPE& first, const VARTYPE& ... arg )
{
  std::vector<OBJTYPE> ans = MakeVector<OBJTYPE>( arg ... );
  ans.insert( ans.begin(), OBJTYPE( first ) );
  return ans;
}

/**
 * @brief Starting point for the recursive call for the variadic interfaces.
 */
template<typename OBJTYPE, typename ARGTYPE>
std::vector<OBJTYPE>
MakeVector( const ARGTYPE& first )
{
  return std::vector<OBJTYPE>( 1, OBJTYPE(first) );
}

/**
 * @brief removing element in a vector is element is equivalent to x
 */
template<typename OBJTYPE>
void ClearValue( std::vector<OBJTYPE>& vec, const OBJTYPE& x )
{
  vec.erase( std::remove( vec.begin(), vec.end(), x ), vec.end() );
}

/**
 * @brief Checking if any of the elements in the vector is equivalent to x
 */
template<typename OBJTYPE>
bool
FindValue( const std::vector<OBJTYPE>& vec, const OBJTYPE& x )
{
  return std::find( vec.begin(), vec.end(), x ) != vec.end();
}

/**
 * @brief Checking if any of the elements in the vector is equivalent to x
 */
template<typename OBJTYPE>
bool
FindValue( std::vector<OBJTYPE>& vec, const OBJTYPE& x )
{
  return std::find( vec.begin(), vec.end(), x ) != vec.end();
}

}/* usr */

#endif/* end of include guard: USERUTILS_COMMON_VECTORUTILS */
