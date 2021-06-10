/**
 * @file    roofit_kstest.cc
 * @brief   Testing functionalities of Proper implementation of KS test using RooFit objects
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#else
#include "UserUtils/MathUtils/RooFitExt.hpp"
#endif

#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "TRandom3.h"

#include <iostream>
using namespace std;

int
main( int argc, char* argv[] )
{
  RooRealVar x( "x", "x", -10, 10 );
  RooRealVar w( "w", "w", -10, 10 );// Just for setting up range? WTF?
  RooDataSet unweighted( "ud", "ud", RooArgSet( x ) );
  RooDataSet weighted( "wd", "wd", RooArgSet( x, w ), RooFit::WeightVar( "w" ) );

  TRandom3 myrand;
  const int N = 2e5;
  for( int i = 0; i < N ; ++i ){
    x = myrand.Gaus( 0, 3 );
    unweighted.add( RooArgSet( x ) );
    x = myrand.Gaus( 0, 3 );
    weighted.add( RooArgSet( x ), myrand.Gaus( 1, 0.1 ) );
  }

  RooGaussian g( "g", "g", x, RooFit::RooConst( 0 ), RooFit::RooConst( 3 ) );

  // cout << usr::KSDistance( weighted, g, x ) << " "
  //      << usr::KSProb( weighted, g, x ) << endl;
  // cout << usr::KSDistance( unweighted, g, x ) << " "
  //      << usr::KSProb( unweighted, g, x ) << endl;
  // cout << usr::KSDistance( weighted, unweighted, x ) << " "
  //      << usr::KSProb( unweighted, weighted, x, true ) << " "
  //      << usr::KSProb( unweighted, weighted, x, false ) << " "
  //      << endl;

  x.setRange("lower",-10,-3);
  x.setRange("upper",3,10);
  x.setRange("center",-3,3);


  cout
      << usr::KSDistance( weighted, g, x, RooFit::CutRange("center") ) << " "
      << usr::KSDistance( weighted, g, x, RooFit::CutRange("lower") ) << " "
      << usr::KSDistance( weighted, g, x, RooFit::CutRange("upper") ) << " "
      << usr::KSDistance( unweighted, g, x, RooFit::CutRange("lower"), RooFit::CutRange("upper") )
  << endl;


  return 0;
}
