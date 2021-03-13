/**
 * @file    StatisticsUtil.cc
 * @brief   Implementation of common statistical number and functions.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/GSLUtil.hpp"
#else
#include "UserUtils/MathUtils/GSLUtil.hpp"
#endif

#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

#include <cassert>
#include <iostream>

namespace usr {

namespace stat {

/**
 * @brief CDF of the normal distribution.
 *
 * Strictly, [GSL's error function](https://www.gnu.org/software/gsl/manual/html_node/Error-Function.html#Error-Function) is defined according to the standard definition of the error
 * function
 *
 * \f[
 *      \mathrm{erf}(x) = \frac{2}{\sqrt(\pi)} \int_0^x dt e^{-t^2}
 * \f]
 *
 * Conversion to CDF of a normal distribution would be:
 *
 * \f[
 *  \mathrm{CDF}(x)
 *      = \frac{1}{2} + \frac{1}{2}\mathrm{erf}\left(\frac{x}{\sqrt{2}}\right)
 * \f]
 */
double
NormalCDF( const double x )
{
  return 0.5 * ( 1 + gsl_sf_erf( x/sqrt( 2. ) ) );
}


/**
 * @brief Translating the informal but common phrase of "x sigmas" into
 *        actual confidence levels.
 */
double
GetConfidenceLevel( const double sigmainterval )
{
  return NormalCDF( sigmainterval ) - NormalCDF( -sigmainterval );
}

/**
 * @brief Translatign the more formal term of "a confidence level of x", into a
 *        more immediately understandable "x sigma interval"
 */
double
GetSigmaInterval( const double confidencelevel )
{
  // Trying to get a larger than 10 sigma is impractical
  static const double maxinterval = 10;

  gsl_root_fsolver* solver = gsl_root_fsolver_alloc( gsl_root_fsolver_brent );

  gsl_function func;
  func.function = []( double x, void* y ){
                    return GetConfidenceLevel( x ) - *( (const double*)y );
                  };
  func.params = const_cast<double*>( &confidencelevel );
  gsl_root_fsolver_set( solver, &func, 0, maxinterval );
  gsl::IterateSolver( solver );

  const double ans = gsl_root_fsolver_root( solver );
  gsl_root_fsolver_free( solver );

  return ans;
}

extern const double onesigma_level = GetConfidenceLevel( 1 );

extern const double twosigma_level = GetConfidenceLevel( 2 );


/**
 * @brief getting the difference in NLL required for a given sigma interval.
 * @details Essentially calculating the difference in NLL from the normal
 * distribution: for a sigma interval of \f$x\f$, the return value would be
 *  \f$x^{2}/2\f$
 */
double
DeltaNLLFromSigma( const double sigma )
{
  return sigma*sigma/2;
}

/**
 * @brief getting the difference in NLL required for a certain confidence level.
 * @details First convert the confidence interval into a sigma interval, then
 * returning the answer.
 */
double
DeltaNLLFromConfidence( const double confidence )
{
  return DeltaNLLFromSigma( GetSigmaInterval( confidence ) );
}

/**
 * @brief Generic routine for calculating minos uncertainties for 1D functions.
 *
 * Given a 1D gsl function acting as the @NLL function for some parameter, the
 * MinosError will use the guess and min/max brackets to return the central
 * value (where the @NLL function was minimized) and the upper and lower
 * uncertainty edges (where the @NLL function deviates from the minimum value
 * by some number K). The number K is determined by the working "confidence
 * level" of the function.
 */
int
MinosError(
  gsl_function* nllfunction,
  double&       guess,
  double&       min,
  double&       max,
  const double  confidencelevel
  )
{
  // Getting central value via minimization
  gsl_min_fminimizer* solver = gsl_min_fminimizer_alloc( gsl_min_fminimizer_brent );
  gsl_min_fminimizer_set( solver, nllfunction, guess, min, max );

  usr::gsl::IterateSolver( solver );

  // saving central value to initial guess variable
  guess = gsl_min_fminimizer_x_minimum( solver );

  gsl_min_fminimizer_free( solver );

  // Getting error via solving nll(x) = nll(x_min)+ sigmainterval^2 / 2
  double deltanll  = DeltaNLLFromConfidence( confidencelevel );
  double nllmin    = nllfunction->function( guess, nllfunction->params );
  double nlltarget = nllmin + deltanll;

  double max_temp = usr::gsl::Solve1D( nllfunction, nlltarget, guess, max );
  double min_temp = usr::gsl::Solve1D( nllfunction, nlltarget, min, guess );

  max = std::max( max_temp, min_temp );
  min = std::min( max_temp, min_temp );

  return 0;
}

/**
 * @brief Generic routine for calculating the minos uncertainty of a N
 *        parameter system.
 *
 * Given an N dimensional parameter system governed by a @NLL function
 * `nllfunction`, the MinosError routine will attempt to calculate the
 * minos uncertainty of the extended parameter defined by the `varfunction`.
 *
 * The `initguess` is mandatory for initializing the parameter guessing, while
 * the `central`, `min`, `max` parameters are no longer used for initializing
 * or bracketting, but rather for storing the return result. The use can
 * provide the parameter inputs for a higher or lower guess if needed, otherwise
 * the guess is simply stepping the minimum value by a value of 0.01.
 */
int
MinosError(
  usr::gsl::gsl_multifunc* nllfunction,
  usr::gsl::gsl_multifunc* varfunction,
  gsl_vector*              initguess,
  double&                  central,
  double&                  min,
  double&                  max,
  const double             confidencelevel,
  gsl_vector*              upperguess,
  gsl_vector*              lowerguess
  )
{
  static const double defaultstep = 0.01;

  // Step 1: Finding the minimum value
  const size_t dim                = nllfunction->n;
  gsl_multimin_fminimizer* solver = gsl_multimin_fminimizer_alloc( gsl_multimin_fminimizer_nmsimplex2, dim );
  gsl_vector* step                = gsl_vector_alloc( dim );
  gsl_vector* bestinput           = gsl_vector_alloc( dim );

  for( size_t i = 0; i < dim; ++i ){
    double stepsize = defaultstep;// Default step size
    const double x  = gsl_vector_get( initguess, i );
    if( upperguess ){
      const double y = gsl_vector_get( upperguess, i );
      stepsize = std::min( stepsize, fabs( y-x )*defaultstep );
    } else if( lowerguess ){
      const double y = gsl_vector_get( upperguess, i );
      stepsize = std::min( stepsize, fabs( y-x )*defaultstep );
    }
    gsl_vector_set( step, i, stepsize );
  }

  gsl_multimin_fminimizer_set( solver, nllfunction, initguess, step );
  usr::gsl::IterateSolver( solver );
  gsl_vector_memcpy( bestinput, solver->x );
  central = varfunction->f( bestinput, varfunction->params );
  double nllmin = nllfunction->f( bestinput, nllfunction->params );

  gsl_vector_free( step );
  gsl_multimin_fminimizer_free( solver );

  // Step 2: Setting up the function for root finding
  // The original variables take up the first 'dim' variables in the solver
  // The final variable is the Lagrange multiplier lambda
  struct minosparam
  {
    usr::gsl::gsl_multifunc* nllfunc;
    usr::gsl::gsl_multifunc* varfunc;
    double                   yval;
  };
  gsl_multiroot_fsolver* errsolver
    = gsl_multiroot_fsolver_alloc( gsl_multiroot_fsolver_hybrids, dim + 1 );

  gsl_multiroot_function minosfunc;
  minosfunc.f
    = []( const gsl_vector* input, void* p, gsl_vector* output ){
        // Getting the original parameters (without the lagrange multipler)
        gsl_vector* reduced = gsl_vector_alloc( input->size - 1 );

        for( size_t i = 0; i < reduced->size; ++i ){
          gsl_vector_set( reduced, i, gsl_vector_get( input, i ) );
        }

        // The lagrange multiplier is the last variable
        double lambda = gsl_vector_get( input, input->size - 1 );

        const minosparam* param = (minosparam*)p;

        // First dim outputs are the partial derivative equations
        for( size_t i = 0; i < reduced->size; ++i ){
          double nlldiv = usr::gsl::partial_deriv( param->nllfunc, reduced, i );
          double vardiv = usr::gsl::partial_deriv( param->varfunc, reduced, i );
          gsl_vector_set( output, i, nlldiv - lambda*vardiv );
        }

        // Final is the NNL value equation
        double nllval = param->nllfunc->f( reduced, param->nllfunc->params );
        gsl_vector_set( output, input->size-1, nllval - param->yval );

        return (int)GSL_SUCCESS;
      };

  minosfunc.n = dim+1;
  minosparam tmp = {
    nllfunction,
    varfunction,
    nllmin + DeltaNLLFromConfidence( confidencelevel )
  };
  minosfunc.params = &tmp;

  // Step 3. Begin the root finding from initial guess
  gsl_vector* errinitguess = gsl_vector_alloc( dim+1 );
  gsl_vector* finalinput   = gsl_vector_alloc( dim );

  for( size_t i = 0; i < dim; ++i ){
    const double setval = upperguess ? gsl_vector_get( upperguess, i ) :
                          gsl_vector_get( bestinput, i )+defaultstep;
    gsl_vector_set( errinitguess, i,  setval );
  }

  gsl_vector_set( errinitguess, dim, 1 );// Initial guess for lambda is 1
  gsl_multiroot_fsolver_set( errsolver, &minosfunc, errinitguess );
  usr::gsl::IterateSolver( errsolver );

  for( size_t i = 0; i < finalinput->size; ++i ){
    gsl_vector_set( finalinput, i, gsl_vector_get( errsolver->x, i ) );
  }

  const double boundary1 = varfunction->f( finalinput, varfunction->params );

  //  Step 3.5 Getting the second root
  for( size_t i = 0; i < dim; ++i ){
    const double best   = gsl_vector_get( bestinput, i );
    const double diff   = gsl_vector_get( finalinput, i ) - best;
    const double setval = lowerguess ? gsl_vector_get( lowerguess, i ) :
                          ( best - diff * 0.01 );
    gsl_vector_set( errinitguess, i, setval );
  }

  // New lambda is the negative of the first lambda
  gsl_vector_set( errinitguess, dim, -gsl_vector_get( errinitguess, dim) );
  gsl_multiroot_fsolver_set( errsolver, &minosfunc, errinitguess );
  usr::gsl::IterateSolver( errsolver );

  for( size_t i = 0; i < finalinput->size; ++i ){
    gsl_vector_set( finalinput, i, gsl_vector_get( errsolver->x, i ) );
  }

  const double boundary2 = varfunction->f( finalinput, varfunction->params );

  // Cleaning up - Storing outputs
  min = std::min( boundary1, boundary2 );
  max = std::max( boundary1, boundary2 );

  gsl_vector_free( errinitguess );
  gsl_vector_free( finalinput   );
  gsl_multiroot_fsolver_free( errsolver );

  return 0;
}


/**
 * @brief GSL compatible function for the @NLL of a Gaussian function.
 * @param x      the estimated parameter value.
 * @param params params[0]=mean, params[1]=sigma
 */
double
GaussianNLL( double x, void* params )
{
  double mean  = ( (double*)params )[0];
  double sigma = ( (double*)params )[1];
  return ( x-mean )*( x-mean )/( 2*sigma*sigma );
}

/**
 * @brief GSL compatible function for the @NLL of a binomial measurement.
 * @param x      the estimated parameter value (\f$\epsilon\f$)
 * @param params param[0] = passed, param[1]=total
 */
double
BinomialNLL( double x, void* params )
{
  double passed = ( (double*)params )[0];
  double total  = ( (double*)params )[1];
  return -passed* log( x ) - ( total-passed ) * log( 1 - x );
}

/**
 * @brief GSL compatible function for the @NLL of a poisson measurement.
 * @param x      the estimated parameter value (\f$n\f$)
 * @param params param[0] = number of observed events.
 */
double
PoissonNLL( double x, void* params )
{
  double obs = ( (double*)params )[0];
  return x - obs*log( x );
}

}/* stat */

}/* usr */
