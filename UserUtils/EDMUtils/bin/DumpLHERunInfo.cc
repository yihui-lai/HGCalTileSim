/**
 * @file    DumpLHERunInfo.cc
 * @ingroup EDMUtils
 * @brief   Dumping contents of the LHE production information on screen.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @details
 * This is meant to be used with a single EDM file:
 * ```
 * DumpLHERunInfo file.root
 * ```
 * The LHE production information of the first run would be printed out.
 * Notice that the "tag" (type of production information) would be printed
 * before each line, so you can combine with grep commands to get specific
 * parts of the production information.
 */

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "TFile.h"

#include <iostream>
using namespace std;

int
main( int argc, char* argv[] )
{
  fwlite::Run run( TFile::Open( argv[1] ) );
  fwlite::Handle<LHERunInfoProduct> lhehandle;

  lhehandle.getByLabel( run, "externalLHEProducer" );

  for( auto iter = lhehandle.ref().headers_begin();
       iter != lhehandle.ref().headers_end(); iter++ ){
    const auto& lines = iter->lines();

    for( const auto line : lines ){
      std::cout << "[" << iter->tag()  << "] " << line << flush;
    }
  }
}
