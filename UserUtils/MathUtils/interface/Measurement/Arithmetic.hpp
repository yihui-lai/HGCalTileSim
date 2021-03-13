/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Defining functions for measurement arithmetics
 */
#ifndef USERUTILS_MATHUTILS_MEASUREMENT_ARITHMETIC_HPP
#define USERUTILS_MATHUTILS_MEASUREMENT_ARITHMETIC_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#include "UserUtils/MathUtils/interface/StatisticsUtil.hpp"
#else
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#include "UserUtils/MathUtils/Measurement/Measurement.hpp"
#include "UserUtils/MathUtils/StatisticsUtil.hpp"
#endif

#include <vector>

namespace usr {

/**
 * @addtogroup StatUtils
 * @{
 */

/*-----------------------------------------------------------------------------
 *  Container interface for  with MinosError functions
   --------------------------------------------------------------------------*/
extern Measurement MakeMinos(
  gsl_function* nll,
  const double  initguess,
  const double  min,
  const double  max,
  const double  confidencelevel = usr::stat::onesigma_level
  );

extern Measurement MakeMinos(
  usr::gsl::gsl_multifunc* nll,
  usr::gsl::gsl_multifunc* varfunction,
  gsl_vector*              initguess,
  const double             confidencelevel = usr::stat::onesigma_level,
  gsl_vector*              upperguess      = nullptr,
  gsl_vector*              lowerguess      = nullptr
  );

/*-----------------------------------------------------------------------------
 *  Default approximation for NLL function for a given measurement.
   --------------------------------------------------------------------------*/
extern double LinearVarianceNLL( double, const Measurement& );

/*-----------------------------------------------------------------------------
 *  Simple calculation of un-correlated function.
   --------------------------------------------------------------------------*/
extern Measurement SumUncorrelated(
  const std::vector<Measurement>& paramlist,
  const double confidencelevel = usr::stat::onesigma_level,
  double ( *nll )( double, const Measurement& ) = & LinearVarianceNLL
  );

extern Measurement ProdUncorrelated(
  const std::vector<Measurement>& paramlist,
  const double confidencelevel = usr::stat::onesigma_level,
  double ( *nll )( double, const Measurement& ) = & LinearVarianceNLL
  );

/**
 * @brief variadic interface to allow for the indefinite input Sum(x,y,z,....)
 *
 * It is worth noting that the output of (x+y+z) might not be the same as the
 * output of Sum(x,y,z) (operator+ is a two input call to Sum(a,b)). But for
 * all intents and purposes, the errors in uncertainties could be ignore for
 * the most part.
 */
template<typename ... Ts>
Measurement
Sum( const Measurement& x,  Ts ... args )
{
  return SumUncorrelated( MakeVector<Measurement>( x, args ... ) );
}

/**
 * @brief variadic interface to allow for the indefinte input Prod(x,y,z,....)
 *
 * It is worth noting that the output of (x*y*z) might not be the same as the
 * output of Prod(x,y,z) (operator* is a two input call to Prod(a,b)). But for
 * all intents and purposes, the errors in uncertainties could be ignore for
 * the most part.
 */
template<typename ... Ts>
Measurement
Prod( const Measurement& x, Ts ... args )
{
  return ProdUncorrelated( MakeVector<Measurement>( x, args ... ) );
}

/** @} */

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_PARAMETERARITHMETIC_HPP */
