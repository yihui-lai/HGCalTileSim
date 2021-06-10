/**
 * @file    Maths.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Simple concrete type maths functions
 */
#ifndef USERUTILS_COMMON_MATHS_HPP
#define USERUTILS_COMMON_MATHS_HPP

#include <cstdint>
#include <vector>

namespace usr {

/**
 * @defgroup STLMaths Math with concrete types.
 * @brief    Calculations with concrete types that are not available in STL
 * @ingroup  Common
 * @{
 */

/**
 * @brief the [sign function](https://en.wikipedia.org/wiki/Sign_function) for
 *        any 0 comparable type.
 *
 * For any type that can be constructed and compared with 0, returns 1 if
 * greater than 0, -1 if smaller than 0, 0 if exactly zero.
 */
template<typename TYPE>
inline int
sgn( TYPE val )
{
  return ( TYPE( 0 ) < val ) - ( val < TYPE( 0 ) );
}

/*-----------------------------------------------------------------------------
 *  exponential and related functions
   --------------------------------------------------------------------------*/
extern int    ReduceToMant( double& x );
extern int    GetExponent( double x );
extern double IntPower( double base, int exp );

/*-----------------------------------------------------------------------------
 *  Rounding functions
   --------------------------------------------------------------------------*/
extern unsigned RoundUpToP2( unsigned x );

/*-----------------------------------------------------------------------------
 *  Hash value calculation
   --------------------------------------------------------------------------*/
extern uint64_t HashValue64( double x );
extern uint32_t HashValue32( double x );
extern uint64_t Hash64Join( uint64_t, uint64_t );
extern uint32_t Hash32Join( uint32_t, uint32_t );

/*-----------------------------------------------------------------------------
 *  Statistics with vector<double>s
   --------------------------------------------------------------------------*/
extern double Mean( const std::vector<double>& );
extern double StdDev( const std::vector<double>& );

/*-----------------------------------------------------------------------------
 *  More rounding functions
   --------------------------------------------------------------------------*/
extern double RoundTo( const double x, const double target );
extern double RoundUp( const double x, const double target );
extern double RoundDown( const double x, const double target );

/** @} */

}/* usr */

#endif