/**
 * @file    StatisticsUtil.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Calculating statistics routines using GSL objects.
 */

#ifndef USERUTILS_MATHUTILS_STATISTICUTIL_HPP
#define USERUTILS_MATHUTILS_STATISTICUTIL_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/GSLUtil.hpp"
#else
#include "UserUtils/MathUtils/GSLUtil.hpp"
#endif

namespace usr {

namespace stat {

/**
 * @addtogroup StatUtils
 * @{
 */

/*-----------------------------------------------------------------------------
 *  Sigma interval -- confidence level conversion.
   --------------------------------------------------------------------------*/
extern double NormalCDF( const double x );
extern double GetConfidenceLevel( const double sigmainterval );
extern double GetSigmaInterval( const double confidencelevel );

extern const double onesigma_level;
extern const double twosigma_level;

/*-----------------------------------------------------------------------------
 *  General calculation routines for Minos uncertainties.
   --------------------------------------------------------------------------*/
extern double DeltaNLLFromSigma( const double sigma );
extern double DeltaNLLFromConfidence( const double confidence );

// Single value Minos error evaluation
extern int MinosError(
  gsl_function* nllfunction,
  double&       guess,
  double&       min,
  double&       max,
  const double  confidencelevel = onesigma_level
  );

extern int MinosError(
  usr::gsl::gsl_multifunc* nllfunction,
  usr::gsl::gsl_multifunc* varfunction,
  gsl_vector*              initguess,
  double&                  central,
  double&                  min,
  double&                  max,
  const double             confidencelevel = onesigma_level,
  gsl_vector*              upperguess = nullptr,
  gsl_vector*              lowerguess = nullptr
  );

/*-----------------------------------------------------------------------------
 *  Common single variable NLL functions
   --------------------------------------------------------------------------*/
extern double GaussianNLL( double, void* params );
extern double BinomialNLL( double, void* params );
extern double PoissonNLL( double, void* params );
/** @} */

}/* stat */

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_STATISTICUTIL_HPP */
