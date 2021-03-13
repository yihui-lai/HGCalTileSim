/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Implementing common measurement constructions
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/Measurement/Arithmetic.hpp"
#include "UserUtils/MathUtils/interface/Measurement/CommonDistro.hpp"
#include "UserUtils/MathUtils/interface/StatisticsUtil.hpp"
#else
#include "UserUtils/MathUtils/Measurement/Arithmetic.hpp"
#include "UserUtils/MathUtils/Measurement/CommonDistro.hpp"
#include "UserUtils/MathUtils/StatisticsUtil.hpp"
#endif

#include "TEfficiency.h"

#include <limits>

namespace usr {

namespace Efficiency {

/**
 * @brief simpler interface for creating Efficiency (binomial) measurements
 *         with minos uncertainties.
 * @ingroup StatUtils
 *
 * @details Calls the GSL functions for generating, known to fail for the
 * extreme case of all pass and all fail.
 */
Measurement
Minos(
  const double passed,
  const double total,
  const double confidencelevel )
{
  double params[] = {passed, total};
  gsl_function binomial;
  binomial.function = &usr::stat::BinomialNLL;
  binomial.params   = params;

  const double est    = passed/total;
  const double min    = usr::gsl::mch_epsilon ;

  return MakeMinos(
    &binomial,
    est,
    min, 1-min,
    confidencelevel
    );
}

const bool shortest_interval = true;
const bool central_interval  = false;

/**
 * @brief Creating a efficiency (binomial) measurement with Bayesian
 *        uncertainties.
 * @ingroup StatUtils
 *
 * @details Using TEfficiency methods, see those for detailed documentation.
 */
Measurement
Bayesian(
  double passed,
  double total,
  double confidencelevel,
  bool   confidencemethod,
  double alpha,
  double beta
  )
{
  const double central = ( passed + alpha - 1. )/( total + alpha + beta - 2. );
  const double err_up  = TEfficiency::Bayesian(
    total, passed,
    confidencelevel,
    alpha, beta,
    true,// For upper boundary
    confidencemethod
    ) - central;
  const double err_down = central - TEfficiency::Bayesian(
    total, passed,
    confidencelevel,
    alpha, beta,
    false,// For lower boundary
    confidencemethod
    );
  return Measurement( central, err_up, err_down );
}

}/* Efficiency */

/*----------------------------------------------------------------------------*/

namespace Poisson {

/**
 * @brief simple interface for creating poisson measurements with minos
 *        uncertainties.
 * @ingroup StatUtils
 */
Measurement
Minos( const double obs, const double confidencelevel )
{
  gsl_function poisson;
  poisson.function = &usr::stat::PoissonNLL;
  poisson.params   = const_cast<double*>( &obs );

  const double min = usr::gsl::mch_epsilon ;
  const double max = obs + obs*obs +1 ;

  return MakeMinos(
    &poisson,
    obs,
    min,
    max,
    confidencelevel );
}

} /* Poisson */

}/* usr */
