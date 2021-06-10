/**
 * @file    format_test.cc
 * @brief   unit testing for format related functions.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @details Detailed descriptions
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/STLUtils.hpp"
#else
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/Common/STLUtils.hpp"
#endif

#include <boost/format.hpp>
#include <iostream>

using namespace std;

void
testdecimal( double x )
{
  cout << "Default: " << usr::fmt::base::decimal( x ) << endl
       << usr::separator() << endl;

  for( int i = 0; i < 30; ++i ){
    usr::fout( "Precision: %3d | %40s | %40s | %40s",
      i,
      usr::fmt::base::decimal( x, i ),
      usr::fmt::base::decimal( x, -i ),
      usr::fmt::base::decimal( x, i ).spacestr( " " ).str() );
  }

  cout << usr::separator() << endl;
}

void
testscientific( double x )
{
  for( int i = 0; i < 30; ++i ){
    usr::fout( "Precision %3d | %40s | %40s\n",
      i,
      usr::fmt::base::scientific( x, i ).str(),
      usr::fmt::base::scientific( x, i ).spacestr( "  " ).str() );
  }

  cout << usr::separator() << endl;
}

int
main( int argc, char* argv[] )
{
  cout << usr::separator() << ">>> General decimal point testing" << endl;

  testdecimal( 123456789 );
  testdecimal( 12345679.123456 );

  testscientific( 123456789 );
  testscientific( 12345679.123456 );

  return 0;
}
