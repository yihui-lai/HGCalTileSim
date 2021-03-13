/**
 * @file    Counter.cc
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Implementing Counter class
 */
#include "UserUtils/EDMUtils/interface/Counter.hpp"

namespace usr {

Counter::Counter(){ value = 0; }
Counter::Counter( const double x ){ value = x; }
Counter::~Counter(){}
/**
 * @brief The existence of this function is required for merging values between
 *        runs.
 *
 * Unlike its edm counterpart. We would not care for overflow parsing here (it is
 * highly unlikely that some cache would require e308 to store.)
 */
bool
Counter::mergeProduct( const Counter& nextcounter )
{
  value += nextcounter.value;
  return true;
}

}/* usr */
