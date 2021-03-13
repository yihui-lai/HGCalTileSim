/**
 * @file    GSLUtil.cc
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Implementing GSL operation functions
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/GSLUtil.hpp"
#else
#include "UserUtils/MathUtils/GSLUtil.hpp"
#endif

#include <cassert>
#include <cmath>
#include <gsl/gsl_deriv.h>
#include <iostream>

namespace usr {

namespace gsl {

/*-----------------------------------------------------------------------------
 *  GSL settings constants
   --------------------------------------------------------------------------*/
/**
 * @brief absolute value difference to terminate iterations.
 */
extern const double abs_epsilon     = 1e-8;

/**
 * @brief relative value difference to terminate iterations.
 */
extern const double rel_epsilon     = 1e-6;

/**
 * @brief pseudo machine epsilon to initialize iterations for next to boundary
 *        situations.
 * @details this is the difference between 1 and the next largest double
 *          smaller than 1 (\f$\sim10^{-50}\f$). Enough for most test cases to
 *          converge properly.
 */
extern const double mch_epsilon     = 1-std::nexttoward( 1, 0 );

/**
 * @brief maximum iterations to run the various GSL solvers.
 * @details should be enough to get decent values for typical use cases.
 */
extern const unsigned max_iteration = 1e+4;

/*-----------------------------------------------------------------------------
 *  Solver iterations.
   --------------------------------------------------------------------------*/
/**
 * @brief automaticaly iterating the gsl single root solver.
 *
 * Notice that this iteration would terminate if the absolute uncertainty
 * is smaller than the gsl::abs_epsilon. So this function might not be suitable
 * is the solution is known to be extremely small.
 */
void
IterateSolver( gsl_root_fsolver* solver )
{
  int status;
  size_t iter = 0;
  double lo;
  double hi;

  do {
    iter++;
    status = gsl_root_fsolver_iterate( solver );
    lo     = gsl_root_fsolver_x_lower( solver );
    hi     = gsl_root_fsolver_x_upper( solver );
    status = gsl_root_test_interval( lo, hi, 0, abs_epsilon );
    if( status == GSL_SUCCESS ){
      break;
    }
  } while( status == GSL_CONTINUE && iter < max_iteration );
}

/**
 * @brief automatically iterating the gsl single variable minimizer.
 *
 * Notice that this iteration would terminate if the absolute uncertainty
 * is smaller than the gsl::abs_epsilon. So this function might not be suitable
 * is the solution is known to be extremely small.
 */
void
IterateSolver( gsl_min_fminimizer* solver )
{
  int status;
  size_t iter = 0;
  double lo;
  double hi;

  do {
    iter++;
    status = gsl_min_fminimizer_iterate( solver );
    lo     = gsl_min_fminimizer_x_lower( solver );
    hi     = gsl_min_fminimizer_x_upper( solver );

    status = gsl_min_test_interval( lo, hi, 0, abs_epsilon );

  } while( status == GSL_CONTINUE && iter < max_iteration );

}

/**
 * @brief automatically handling the gsl solver for a multiroot setups
 *
 * The maximum number of iterations is not just gsl::max_iteration, but timed by
 * the dimension of the problem given in the solver. As the residual test only
 * provided by GSL only contains distance tests, the absolute distance epsilon
 * is used for termination determination.
 */
void
IterateSolver( gsl_multiroot_fsolver* solver )
{
  const unsigned dim                = solver->f->size;
  const unsigned this_max_iteration = max_iteration * dim;
  int status;
  size_t iteration = 0;

  do {
    iteration++;
    status = gsl_multiroot_fsolver_iterate( solver );
    // if( status ){ break; }
    status = gsl_multiroot_test_residual( solver->f, abs_epsilon );
  } while( status == GSL_CONTINUE && iteration < this_max_iteration );
}

/**
 * @brief automatically iterating a gsl multi-dimension minimizer.
 *
 * The maximum number of iterations is not just gsl::max_iteration, but timed by
 * the dimension of the problem given in the solver. As the residual test only
 * provided by GSL only contains distance tests, the absolute distance epsilon
 * is used for termination determination.
 */
void
IterateSolver( gsl_multimin_fminimizer* solver )
{
  int status;
  unsigned iter = 0;
  double size;// size of stepping vector

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate( solver );

    if( status ){ break; }

    size   = gsl_multimin_fminimizer_size( solver );
    status = gsl_multimin_test_size( size, abs_epsilon );

  } while( status == GSL_CONTINUE && iter < max_iteration * solver->x->size );
}

/**
 * @brief solving a 1 dimension function \f$f(x)=y\f$, user must provide bracket
 */
double
Solve1D(
  gsl_function* func,
  double        yval,
  double        xmin,
  double        xmax
  )
{
  struct argparam
  {
    gsl_function* original_func;
    double        yval;
  };
  argparam newparams = { func, yval };

  gsl_function augfunc;
  augfunc.function = []( double x, void* param ){
                       argparam* p = (argparam*)param;
                       return p->original_func->function( x, p->original_func->params ) - p->yval;
                     };
  augfunc.params = &newparams;

  gsl_root_fsolver* solver = gsl_root_fsolver_alloc( gsl_root_fsolver_brent );
  gsl_root_fsolver_set( solver, &augfunc, xmin, xmax );

  IterateSolver( solver );

  double ans = gsl_root_fsolver_root( solver );

  gsl_root_fsolver_free( solver );

  return ans;
}

/**
 * @brief handy parameter-less function for later uses.
 */
double
sum( const gsl_vector* x, void* )
{
  double ans = 0;

  for( size_t i = 0; i < x->size; ++i ){
    ans += gsl_vector_get( x, i );
  }

  return ans;
}

/**
 * @brief handy parameter-less function for later use.
 */
double
product( const gsl_vector* x, void* )
{
  double ans = 1;

  for( size_t i = 0; i < x->size; ++i ){
    ans *= gsl_vector_get( x, i );
  }

  return ans;
}

/**
 * @brief partial derivative function for a multi-variable function
 *
 * Given the function \f$f(\vec{x})\f$, a input point \f$\vec{x}_0\f$, and a
 * variable index $i$, this function calculates the results
 * \f$\frac{\partial f(\vec{x}_0)}{\partial x_i}\f$.
 * This is using the vanilla GSL interface where the user must provide a
 * stepsize and references to store the results.
 */
int
partial_deriv(
  gsl_multifunc* function,
  gsl_vector*    x,
  size_t         varidx,
  double         h,
  double&        result,
  double&        error
  )
{
  struct projected_param
  {
    gsl_multifunc* original_func;
    gsl_vector*    original_x;
    size_t         varidx;
  };

  gsl_function f;

  f.function =
    []( double x, void* p )
    {
      projected_param* param = (projected_param*)p;

      gsl_vector* newvec = gsl_vector_alloc( param->original_x->size );

      for( size_t idx = 0; idx < param->original_x->size; ++idx ){
        const double newx = idx==param->varidx ?  x:
                            gsl_vector_get( param->original_x, idx );
        gsl_vector_set( newvec, idx, newx );
      }

      double ans = param->original_func->f(
        newvec, param->original_func->params
      );

      gsl_vector_free( newvec );

      return ans;
    };

  projected_param p = { function, x, varidx };

  f.params = &p;

  double projx = gsl_vector_get( x, varidx );
  int retans   = gsl_deriv_central( &f, projx, h, &result, &error );
  return retans;
}

/**
 * @brief simpler interface for calculating a derivative at a given point.
 *
 * The stepsize is automatically determined such that the the error is smaller
 * than the abs_epsilon and rel_epsilon values. The step size starts at
 * \f$10^{-4}\f$ and reduces by a factor of 2 with a minimum cap at
 * \f$10^{-8}\f$
 */
double
deriv( gsl_function* F, double x )
{
  double stepsize = sqrt(abs_epsilon);
  double result = 0, error = 1;

  do {
    gsl_deriv_central( F, x, stepsize, &result, &error );
    stepsize /= 2;
  } while(
    fabs( error ) > std::max( abs_epsilon, result*rel_epsilon ) &&
    stepsize > abs_epsilon
    );

  return result;
}

/**
 * @brief simpler interface for partial derivatives.
 *
 * same scheme as that used in deriv()
 */
double
partial_deriv(
  gsl_multifunc* function,
  gsl_vector*    x,
  size_t         varidx
  )
{
  double stepsize = sqrt(abs_epsilon);
  double result = 0, error = 1;

  do {
    partial_deriv( function, x, varidx, stepsize, result, error );
    stepsize /= 2;
  } while(
    fabs( error ) > std::max( abs_epsilon, result*rel_epsilon ) &&
    stepsize > abs_epsilon
    );

  return result;
}

};/* namespace gsl  */

};/* namespace usr */
