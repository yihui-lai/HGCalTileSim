/**
 * @file    GSLUtil.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Listing GSL operation functions and variables.
 */

#ifndef USERUTILS_MATHUTILS_GSLSETTINGS_HPP
#define USERUTILS_MATHUTILS_GSLSETTINGS_HPP

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_min.h>

namespace usr {

namespace gsl {
/**
 * @defgroup GSLUtils GSLUtils
 * @brief   Default settings and routines for using GSL objects.
 * @ingroup MathUtils
 * @details
 * Simplifying commonly used GSL routines used in various calculation routines
 * to single functions that should converge in for most of the use cases.
 * Extreme cases would always need specific processing, but this collections
 * of functions is still useful for development.
 * @{
 */

/*-----------------------------------------------------------------------------
 *  GSL setting constants
   --------------------------------------------------------------------------*/
extern const double abs_epsilon;
extern const double rel_epsilon;
extern const double mch_epsilon;
extern const unsigned max_iteration;

/*-----------------------------------------------------------------------------
 *  Default GSL iterators
   --------------------------------------------------------------------------*/
extern void IterateSolver( gsl_multiroot_fsolver* solver );
extern void IterateSolver( gsl_root_fsolver* solver );
extern void IterateSolver( gsl_min_fminimizer* solver );
extern void IterateSolver( gsl_multimin_fminimizer* solver );

/*-----------------------------------------------------------------------------
 *  Auto-iteration for solving one-dimensional function at f(x)=a
   --------------------------------------------------------------------------*/
extern double Solve1D(
  gsl_function* function,
  double        yval,
  double        xmin,
  double        ymax
);

/**
 * @brief defining a new name for the \f$R^{n}\rightarrow R\f$ function.
 * @details basically the function used for the gsl minimizers:
 * - double (* f)( const gsl_vector* x, void* params)
 * - size_t n
 * - void* params
 */
typedef gsl_multimin_function gsl_multifunc;

/*-----------------------------------------------------------------------------
 *  Handy parameter-less mult variable functions.
   --------------------------------------------------------------------------*/
double sum( const gsl_vector* x, void* params);
double product( const gsl_vector* x, void* params);

/*-----------------------------------------------------------------------------
 *  Partial derivative function
   --------------------------------------------------------------------------*/
extern int partial_deriv(
  gsl_multifunc* function,
  gsl_vector*    x,
  size_t         varidx,
  double         h,
  double&        result,
  double&        error
  );

extern double deriv(
  gsl_function* function,
  double        x
);

extern double partial_deriv(
  gsl_multifunc* function,
  gsl_vector*    x,
  size_t         varidx
);

/** @} */

};/* gsl */

}/* usr */

#endif/* end of include guard: USERUTILS_MATHUTILS_GSLSETTINGS_HPP */
