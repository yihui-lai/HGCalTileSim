/*******************************************************************************
*
*  Filename    : computing_test.cc
*  Description : Unit testing functions for functions defined in Compute.hpp
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*******************************************************************************/
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement/Format.hpp"
#include "UserUtils/MathUtils/interface/StatisticsUtil.hpp"
#else
#include "UserUtils/Common/STLUtils.hpp"
#include "UserUtils/MathUtils/Measurement/Format.hpp"
#include "UserUtils/MathUtils/StatisticsUtil.hpp"
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace usr;
int
main( int argc, char* argv[] )
{
  {// Constants checking
    cout << separator() << endl
         << "Confidence level from sigma " << endl
         << fmt::decimal( stat::GetConfidenceLevel( 1 ), 8 )<< endl
         << fmt::decimal( stat::GetConfidenceLevel( 2 ), 8 )<< endl;

    cout << separator() << endl
         << "Sigma value from confidence level" << endl
         << fmt::decimal( stat::GetSigmaInterval( 0.68 ), 8 ) << endl
         << fmt::decimal( stat::GetSigmaInterval( 0.95 ), 8 ) << endl
         << fmt::decimal( stat::GetSigmaInterval( stat::onesigma_level ), 8 ) << endl
         << fmt::decimal( stat::GetSigmaInterval( stat::twosigma_level ), 8 ) << endl;
  }

  {// Partial derivative test
    cout << separator() << endl
         << "Partial Derivative test" << endl;

    gsl::gsl_multifunc func;
    func.f = []( const gsl_vector* x, void* ){
               double ans = 0;

               for( size_t i = 0; i < x->size; ++i ){
                 ans += gsl_vector_get( x, i ) * gsl_vector_get( x, i );
               }

               return ans;
             };
    func.n      = 3;
    func.params = nullptr;

    gsl_vector* x = gsl_vector_alloc( 3 );
    gsl_vector_set( x, 0, 1 );
    gsl_vector_set( x, 1, 2 );
    gsl_vector_set( x, 2, 3 );
    cout << func.f( x,  nullptr ) << endl;

    for( double i = 0; i < 5; i += 0.001 ){
      gsl_vector_set( x, 0, i );
      const double div = usr::gsl::partial_deriv( &func, x, 0 );
      const double ext = 2 * i;
      cout <<  fmt::decimal( div-ext, 8 ) << endl;
    }
  }

  {// 1D minor error test
    cout << separator() << endl
         << "1D Minos Error testing" << endl;
    gsl_function gaussnll;
    gaussnll.function = []( double x, void* p ){
                          double* param = (double*)p;
                          double mean   = param[0];
                          double sig    = param[1];
                          return ( x-mean )*( x-mean )/( 2.*sig*sig );
                        };
    {
      double test[] = { 1, 2};
      double mean = 0, errup = 10, errlo = -10;

      gaussnll.params = test;
      usr::stat::MinosError( &gaussnll, mean, errlo, errup );
      cout << mean << " "<< errup  << " " << errlo << endl;
      cout << gaussnll.function( mean, test ) << endl;
      cout << gaussnll.function( errup, test ) << endl;
      cout << gaussnll.function( errlo, test ) << endl;
    }
  }

  {// 1D minos error test
    cout << separator() << endl
         <<"2D minos error test" << endl;
    usr::gsl::gsl_multifunc twogaussnll;
    twogaussnll.f = []( const gsl_vector* x, void* param ){
                      double* p    = (double*)param;
                      double mean1 = p[0];
                      double sig1  = p[1];
                      double mean2 = p[2];
                      double sig2  = p[3];
                      double x1    = gsl_vector_get( x, 0 );
                      double x2    = gsl_vector_get( x, 1 );
                      return ( x1-mean1 )*( x1-mean1 )/( 2.*sig1*sig1 )
                             +( x2-mean2 )*( x2-mean2 )/( 2.*sig2*sig2 );
                    };
    twogaussnll.n = 2;
    double test[] = {3, 0.5, 5, 1};
    twogaussnll.params = test;

    usr::gsl::gsl_multifunc twosum;
    twosum.f = []( const gsl_vector* x, void* ){
                 double ans = 0;

                 for( size_t i = 0; i < x->size; ++i ){
                   ans += gsl_vector_get( x, i );
                 }

                 return ans;
               };
    twosum.n = 2;

    gsl_vector* init = gsl_vector_alloc( 2 );
    gsl_vector_set_all( init, 1 );

    double central, up, down;
    usr::stat::MinosError( &twogaussnll, &twosum, init, central, up, down );
    cout << central << "  " << up << " " << down << endl;
  }

  {// Predefined NLL test
    cout << separator() << endl
         <<"Pre-defined  NLL test" << endl;
    double a           = 1;
    double b           = -100;
    double c           = 100;
    double gausparam[] = {5, 2};
    gsl_function test  = {&usr::stat::GaussianNLL, &gausparam};
    usr::stat::MinosError( &test, a, b, c );
    cout << "Gaussian" << endl;
    cout << fmt::decimal( b, 8 ) << " "
         << fmt::decimal( a, 8 ) << " "
         << fmt::decimal( c, 8 ) <<endl;

    {
      const double passed    = 0.01;
      const double total     = 200;
      double binomialparam[] = {passed, total};
      a    = passed/total;
      b    = gsl::mch_epsilon;
      c    = 1-b;
      test = {&usr::stat::BinomialNLL, &binomialparam};
      usr::stat::MinosError( &test, a, b, c );
      cout << "Binomial" << endl;
      cout << fmt::decimal( b, 8 ) << " "
           << fmt::decimal( a, 8 ) << " "
           << fmt::decimal( c, 8 ) <<endl;
    }

    {
      const double passed    = 1;
      const double total     = 2000000;
      double binomialparam[] = {passed, total};
      a    = passed/total;
      b    = gsl::mch_epsilon;
      c    = 1-b;
      test = {&usr::stat::BinomialNLL, &binomialparam};
      usr::stat::MinosError( &test, a, b, c );
      cout << "Binomial" << endl;
      cout << fmt::decimal( b, 8 ) << " "
           << fmt::decimal( a, 8 ) << " "
           << fmt::decimal( c, 8 ) <<endl;
    }
    {
      double poissonparam = 0.05;
      a    = poissonparam;
      b    = 1-std::nexttoward( 1, 0 );
      c    = a + a*a +  1;
      test = {&usr::stat::PoissonNLL, &poissonparam};
      usr::stat::MinosError( &test, a, b, c );
      cout << "Poisson" << endl;
      cout << fmt::decimal( b, 8 ) << " "
           << fmt::decimal( a, 8 ) << " "
           << fmt::decimal( c, 8 ) <<endl;
    }
    {
      double poissonparam = 1e6;
      a    = poissonparam;
      b    = 1-std::nexttoward( 1, 0 );
      c    = a + a*a +  1;
      test = {&usr::stat::PoissonNLL, &poissonparam};
      usr::stat::MinosError( &test, a, b, c );
      cout << "Poisson" << endl;
      cout << fmt::decimal( b, 8 ) << " "
           << fmt::decimal( a, 8 ) << " "
           << fmt::decimal( c, 8 ) <<endl;
    }
  }

  return 0;
}
