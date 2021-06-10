/**
 * @file    Measurement_Arithmetic.cc
 * @brief   Implementation of measurement arithetics methods
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#include "UserUtils/MathUtils/interface/GSLUtil.hpp"
#include "UserUtils/MathUtils/interface/StatisticsUtil.hpp"
#else
#include "UserUtils/MathUtils/Measurement.hpp"
#include "UserUtils/MathUtils/GSLUtil.hpp"
#include "UserUtils/MathUtils/StatisticsUtil.hpp"
#endif

#include <cassert>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <vector>

using namespace std;

namespace usr {

/**
 * @brief interface to the Single variable version of the MinosError() function.
 */
Measurement
MakeMinos(
  gsl_function* nll,
  const double  initguess,
  const double  min,
  const double  max,
  const double  confidencelevel
  )
{
  double central    = initguess;
  double lowerbound = min;
  double upperbound = max;

  usr::stat::MinosError( nll, central, lowerbound, upperbound, confidencelevel );
  return Measurement( central, upperbound-central, central-lowerbound );
}

/**
 * @brief interface to the \f$R^{n}\rightarrow R\f$ version of the MinosError()
 *        function
 */
Measurement
MakeMinos(
  usr::gsl::gsl_multifunc* nll,
  usr::gsl::gsl_multifunc* varfunction,
  gsl_vector*              initguess,
  const double             confidencelevel,
  gsl_vector*              upperguess,
  gsl_vector*              lowerguess
  )
{
  double central;
  double upperbound;
  double lowerbound;

  usr::stat::MinosError(
    nll,
    varfunction,
    initguess,
    central,
    lowerbound,
    upperbound,
    confidencelevel,
    upperguess,
    lowerguess
    );

  return Measurement( central, upperbound-central, central - lowerbound );
}


static double
AugmentedDeminator( double x, const double up, const double lo )
{
  // Minimum product of both uncertainties.
  // Here we are taking the cut off absolute uncertainty used in the GLS
  // functions.
  static const double minprod = gsl::rel_epsilon* gsl::rel_epsilon;

  const double V = std::max( up * lo, minprod );
  const double A = ( up-lo ) / V;

  double ans = V * ( 1 + x * A );
  if( up > lo && A != 0 && x < ( -lo-1/A )/ 2 ){
    const double s = ( -lo-1/A )/2;
    const double b = A / ( 1+A*s );
    const double B = V * ( 1+A*s ) / exp( b*s );
    ans = B * exp( b * x );
  } else if( lo > up && A != 0 && x > ( up+1/A )/2 ){
    ans = AugmentedDeminator( -x, lo, up );
  }

  return ans;
}

/**
 * @brief Approximate @NLL function for a measurements of
 *        \f$x^{+\sigma_+}_{-\sigma_i}\f$
 *
 * The approximation for an @NLL function is a variation of the one recommended
 * in
 * [R. Barlow's "Asymmetric statistical errors"](https://arxiv.org/abs/physics/0406120),
 * In the form of:
 *
 * \f[
 *  \mathrm{NLL}(x) = \frac{x^{2}}{ 2V(1+Ax)}
 * \f]
 *
 * with \f$V\f$ being the average variance \f$\sigma_+ \sigma_-\f$, and \f$A\f$
 * being the variance asymmetry \f$ (\sigma_{+} - \sigma_{-})/V \f$.
 *
 * A few adjustments are made to ensure the stability of the calculations.
 * - The ratio between the two uncertainties would have a hard limit of 10. The
 *   smaller uncertainty would be inflated if extreme asymmetries are present.
 * - The average variance would have an minimum value of $10^{-16}$ which
 *   should cover most cases where uncertainty calculations are required. Test
 *   have shown that the calculation would not be affected by such limits if
 *   the there are some dominate uncertainties.
 * - To avoid the singularity at \f$x=-1/A\f$, an exponential function would be
 *   used instead for the linear function for the denominator in the @NLL
 *   function when x is past half-way from the central value to the original
 *   singular point, such that the approximation would be extended to an
 *   infinite domain. The exponential function is designed such that the
 *   denominator and its derivative is continuous.
 */
double
LinearVarianceNLL( double x, const Measurement& parm )
{
  static const double maxrelerror = 10.;

  // Getting potentially inflating uncertainties
  const double central  = parm.CentralValue();
  const double rawupper = parm.AbsUpperError();
  const double rawlower = parm.AbsLowerError();
  const double effupper = std::max( rawupper, rawlower / maxrelerror );
  const double efflower = std::max( rawlower, rawupper / maxrelerror );

  const double num = ( x- central ) * ( x - central );
  const double de  = AugmentedDeminator( x-central, effupper, efflower );
  const double ans = 0.5  * num / de;

  return ans;
}

/*-----------------------------------------------------------------------------
 *  Hidden method for making a master NLL function
   --------------------------------------------------------------------------*/
struct uncorrelatedparam
{
  const vector<Measurement>* listptr;
  double                     (* nll)( double, const Measurement& );
};

static double
UncorrelatedNLL( const gsl_vector* x, void* param )
{
  uncorrelatedparam* p = (uncorrelatedparam*)param;
  double ans           = 0;

  for( size_t i = 0; i < p->listptr->size(); ++i ){
    const Measurement& thisp = p->listptr->at( i );
    const double thisx       = gsl_vector_get( x, i );
    ans += p->nll( thisx, thisp );
  }

  return ans;
}

/**
 * @brief given a list of measurements with uncertainties, return the effective
 *        sum of all measurements as if all measurements are uncorrelated.
 *
 * This function also allows the user to define the working confidence level
 * for the calculations and also the approximate @NLL function used for the
 * individual measurements (by default using the LinearVarianceNLL() function).
 *
 * For improved stability of the calculation. An estimate for when the
 * original estimations would give the "1 sigma" uncertainty is given using
 * the symmetric-uncorrelated sum approximation.
 */
Measurement
SumUncorrelated(
  const vector<Measurement>& paramlist,
  const double confidencelevel,
  double ( *nll )( double, const Measurement& )
  )
{
  const unsigned dim = paramlist.size();

  usr::gsl::gsl_multifunc masternll;
  uncorrelatedparam param = { &paramlist, nll };

  masternll.f      = &UncorrelatedNLL;
  masternll.n      = dim;
  masternll.params = &param;

  usr::gsl::gsl_multifunc varfunction;
  varfunction.f      = &usr::gsl::sum;
  varfunction.n      = dim;
  varfunction.params = nullptr;// parameter-less

  gsl_vector* initguess  = gsl_vector_alloc( dim );
  gsl_vector* upperguess = gsl_vector_alloc( dim );
  gsl_vector* lowerguess = gsl_vector_alloc( dim );

  double sqrt_sumsq_up = 0;
  double sqrt_sumsq_lo = 0;

  for( size_t i = 0; i < paramlist.size(); ++i ){
    const Measurement& p = paramlist.at( i );
    gsl_vector_set( initguess, i, p.CentralValue() );
    sqrt_sumsq_up += p.AbsUpperError() * p.AbsUpperError();
    sqrt_sumsq_lo += p.AbsLowerError() * p.AbsLowerError();
  }

  // Using symmetric approximation to get initial guess
  sqrt_sumsq_lo = sqrt( sqrt_sumsq_lo );
  sqrt_sumsq_up = sqrt( sqrt_sumsq_up );

  for( size_t i = 0; i < paramlist.size(); ++i ){
    const Measurement& p = paramlist.at( i );
    const double shift_up = sqrt_sumsq_up / p.AbsUpperError();
    const double shift_lo = sqrt_sumsq_lo / p.AbsLowerError();
    gsl_vector_set( upperguess, i,
      p.CentralValue() + p.AbsUpperError()/shift_up );
    gsl_vector_set( lowerguess, i,
      p.CentralValue() - p.AbsLowerError()/shift_lo );
  }

  // Calling MinosError function for actual computation
  Measurement ans = MakeMinos(
    &masternll,
    &varfunction,
    initguess,
    confidencelevel,
    upperguess,
    lowerguess
    );

  // Releasing gsl vector calls
  gsl_vector_free( initguess );
  gsl_vector_free( upperguess );
  gsl_vector_free( lowerguess );

  return ans;
}

/**
 * @brief given a list of measurements with uncertainties, return the effective
 *        sum of all measurements as if all measurements are uncorrelated.
 *
 * This function also allows the user to define the working confidence level
 * for the calculations and also the approximate @NLL function used for the
 * individual measurements (by default using the LinearVarianceNLL() function).
 */
Measurement
ProdUncorrelated(
  const std::vector<Measurement>& paramlist,
  const double confidencelevel,
  double ( *nll )( double, const Measurement& )
  )
{
  const unsigned dim = paramlist.size();

  // Getting normalized version of list
  double prod = 1.;
  vector<Measurement> normlist;

  for( const auto& p : paramlist ){
    prod *= p.CentralValue();
    normlist.push_back( p.NormParam() );
  }

  usr::gsl::gsl_multifunc masternll;
  uncorrelatedparam param = { &normlist, nll };

  masternll.f      = &UncorrelatedNLL;
  masternll.n      = dim;
  masternll.params = &param;

  usr::gsl::gsl_multifunc varfunction;
  varfunction.f      = &usr::gsl::product;
  varfunction.n      = dim;
  varfunction.params = nullptr;// parameter-less

  gsl_vector* initguess  = gsl_vector_alloc( dim );
  gsl_vector* upperguess = gsl_vector_alloc( dim );
  gsl_vector* lowerguess = gsl_vector_alloc( dim );

  for( size_t i = 0; i < normlist.size(); ++i ){
   const Measurement& p = normlist.at( i );
   gsl_vector_set( initguess,  i, p.CentralValue() );
   // Don't attempt to scale the guess, works better for products
   gsl_vector_set( upperguess, i, p.CentralValue()+p.AbsUpperError() );
   gsl_vector_set( lowerguess, i, p.CentralValue()-p.AbsLowerError() );
  }

  // Calling MinosError function for actual computation
  Measurement ans = MakeMinos(
    &masternll,
    &varfunction,
    initguess,
    confidencelevel,
    upperguess,
    lowerguess
    );

  // Releasing gsl vector calls
  gsl_vector_free( initguess );
  gsl_vector_free( upperguess );
  gsl_vector_free( lowerguess );

  return ans * prod;
}

}/* usr */
