/**
 * @file    Miscellaneous.hpp
 * @brief   Miscellaneous mathematical functions.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#ifndef USERUTILS_MATHUTILS_COMMON_HPP
#define USERUTILS_MATHUTILS_COMMON_HPP

#include "TMatrixD.h"
#include "TMatrixDSym.h"

namespace usr {

/**
 * @defgroup mathmisc Miscellaneous maths functions.
 * @brief A bunch of functions that I find useful but don't know how to classify
 * @ingroup MathUtils
 * @{
 */
extern bool Intersect(
  const double p1x, const double p1y,
  const double p2x, const double p2y,
  const double p3x, const double p3y,
  const double p4x, const double p4y,
  double& intx, double& inty
  );

extern TMatrixD DecompCorvariance( const TMatrixDSym& m );

/** @} */

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_COMMON_HPP */
