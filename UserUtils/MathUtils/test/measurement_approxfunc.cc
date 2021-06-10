/**
 * @file    parameter_approxfunc_test
 * @brief   testing the functionality of likelihood function approximations
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

#include <cstdio>
#include <cassert>
#include <gsl/gsl_deriv.h>

double f( double x, void* p )
{
  return usr::LinearVarianceNLL( x , *((usr::Measurement*)p) ) ;
}

int main(int argc, char* argv[])
{
  usr::Measurement p(1,0.00001,0.00001);
  gsl_function F;
  F.function = &f;
  F.params   = &p;

  double deriv;

  for( double x = -5 ; x < 5 ; x += 0.01 ){
    deriv = usr::gsl::deriv( &F, x );
    printf("%.6lf: %.6lf %.6lf\n", x, f(x,&p), deriv );
    assert( !isnan(deriv) );
  }

  return 0;
}
