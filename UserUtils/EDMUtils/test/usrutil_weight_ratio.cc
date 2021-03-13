/**
 * @file    weight_ratio.cc
 * @brief   Testing the weight ratio information stored in the LHEEventProducts
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "TFile.h"

#include <iostream>
#include <map>
#include <string>
using namespace std;

int
main( int argc, char* argv[] )
{
  fwlite::Event evt( TFile::Open( argv[1] ) );
  fwlite::Handle<LHEEventProduct> lhehandle;
  map<string, double> weightlist ;

  lhehandle.getByLabel( evt, "externalLHEProducer" );

  for( evt.toBegin() ; !evt.atEnd() ; ++evt ){
    weightlist["hepup"] += lhehandle.ref().hepeup().XWGTUP ;
    weightlist["original"] += lhehandle.ref().originalXWGTUP();
    for( const auto& weightinfo : lhehandle.ref().weights() ) {
      weightlist[weightinfo.id] += weightinfo.wgt;
    }
  }

  for( const auto& p : weightlist ){
    printf("%s: %.6lf\n", p.first.c_str(), p.second );
  }
}
