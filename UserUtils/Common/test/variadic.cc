/*******************************************************************************
*
*  Filename    : variadic_test.cc
*  Description : Testing variadic template functions
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*******************************************************************************/
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils.hpp"
#else
#include "UserUtils/Common/STLUtils.hpp"
#endif

#include <vector>
#include <iostream>

using namespace std;

int
main( int argc, char const* argv[] )
{
  vector<int> test = usr::MakeVector<int>(1,2,3,4,5);
  vector<int> test1 = usr::MakeVector<int>(1);
  vector<int> test2 = usr::MakeVector<int>(1,2);
  vector<int> test3 = usr::MakeVector<int>(1,2);

  cout << test.size() << endl;
  for( auto x : test ){
    cout << x << endl;
  }

  cout << usr::fstr( "this is a test" ) << std::endl;
  cout << usr::fstr( "My test %d\n", 1234) << std::endl;

  usr::fout( "this is a test\n" );
  usr::fout( "My test %d\n", 1234) << std::endl;


  return 0;
}
