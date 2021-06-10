/**
 * @file
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Functions for constructing common measurements
 */
#ifndef USERUTILS_MATHUTILS_MEASUREMENT_COMMONDISTRO_HPP
#define USERUTILS_MATHUTILS_MEASUREMENT_COMMONDISTRO_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/Measurement/Measurement.hpp"
#include "UserUtils/MathUtils/interface/StatisticsUtil.hpp"
#else
#include "UserUtils/MathUtils/Measurement/Measurement.hpp"
#include "UserUtils/MathUtils/StatisticsUtil.hpp"
#endif

namespace usr {

/*-----------------------------------------------------------------------------
 *  Efficiency related structures and helper functions
   --------------------------------------------------------------------------*/
namespace Efficiency
{

extern Measurement Minos(
  const double passed,
  const double total,
  const double confidencelevel = usr::stat::onesigma_level
  );

extern const bool shortest_interval;
extern const bool central_interval;

extern Measurement Bayesian(
  double passed,
  double total,
  double confidencelevel =usr::stat::onesigma_level,
  bool   confidencemethod = shortest_interval,
  double alpha = 0.5,
  double beta  = 0.5
  );

};// Efficiency


/*-----------------------------------------------------------------------------
 *  Poisson measurement construction
   --------------------------------------------------------------------------*/
namespace Poisson
{
extern Measurement Minos(
  const double obs ,
  const double confidencelevel = usr::stat::onesigma_level
);

// TODO: Add more methods

};// Poisson

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_EFFICIENCY_HPP */
