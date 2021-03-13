/**
 * @file    Counter.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Defning the Counter class
 */

#ifndef USERUTILS_EDMUTILS_COUNTER_HPP
#define USERUTILS_EDMUTILS_COUNTER_HPP

namespace usr {

/**
 * @brief Class exposed to the EDM ROOT framework to allow for run-level
 *        caching.
 * @ingroup EDMUtils
 *
 * @details
 * This class is essentially the same as the [`edm::MergeableCounter`](http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_7/doc/html/dd/d00/structedm_1_1MergeableCounter.html) class except with doubles
 * instead of integer values.
 */
class Counter
{
public:
  Counter();
  Counter( const double );
  virtual
  ~Counter();

  bool mergeProduct( const Counter& nextCounter );

  // Saving value
  double value;
};

}/* usr  */

#endif/* end of include guard: USERUTILS_EDMUTILS_COUNTER_HPP */
