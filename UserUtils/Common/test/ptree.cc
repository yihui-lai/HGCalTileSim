/*******************************************************************************
*
*  Filename    : ptree_test.cc
*  Description : Unit testing for pt::ptree utility functions
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*******************************************************************************/
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/interface/STLUtils.hpp"
#else
#include "UserUtils/Common/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/STLUtils.hpp"
#endif

#include <iostream>

using namespace std;

/******************************************************************************/

int
main( int argc, char const* argv[] )
{
  cout << usr::separator() << endl
       << "test.json input" << endl;
  usr::pt::ptree tree = usr::FromJsonFile( "test.json" );

  usr::PrintPTree( tree );

  cout << usr::separator() << endl;

  cout << usr::GetSingle<string>( tree, "Static1" ) << endl;

  cout << usr::separator() << endl;

  for( const auto x : usr::GetList<int>( tree, "Instance1.var1" ) ){
    cout << x << endl;
  }

  cout << usr::separator() << endl;

  usr::pt::ptree subtree = usr::GetSubTree( tree, "Instance2" );
  usr::PrintPTree( subtree );

  cout << usr::separator() << endl
       << "name1.json input" << endl;

  usr::pt::ptree tree2 = usr::FromJsonFile( "name1.json" );
  subtree = usr::GetSubTree( tree2, "channel.MuonSignal" );
  cout << usr::CheckQuery( tree2, "channel.MuonSignal" ) << endl;
  cout << usr::CheckQuery( tree2, "channel.MuonSignal.Root Name" ) << endl;
  cout << usr::CheckQuery( tree2, "channel.Garbage" ) << endl;
  cout << usr::CheckQuery( subtree, "Root Name" ) << endl;
  usr::PrintPTree( subtree );

  return 0;
}
