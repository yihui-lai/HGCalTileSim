/*******************************************************************************
*
*  Filename    : parameter_arithmatics.cc
*  Description : Testing the arithmatics of paramters
*  Author      : Yi-Mu "Enoch" Chen [ ensc@hep1.phys.ntu.edu.tw ]
*
*******************************************************************************/
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils.hpp"
#include "UserUtils/MathUtils/interface/Measurement.hpp"
#else
#include "UserUtils/Common/STLUtils.hpp"
#include "UserUtils/MathUtils/Measurement.hpp"
#endif

#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
using namespace usr;

void
AddTest( const vector<Measurement>& list )
{
  std::stringstream ss;

  for( size_t i = 0; i < list.size(); ++i ){
    ss << list[i].CentralValue();
    if( i != list.size() -1 ){
      ss << "+";
    }
  }

  const Measurement sum = SumUncorrelated( list );

  usr::fout( "%40s & %15s & %15s & %15s",
    ss.str(),
    fmt::decimal( sum.CentralValue(),  8 ),
    fmt::decimal( sum.AbsUpperError(), 8 ),
    fmt::decimal( sum.AbsLowerError(), 8 ) );
}


void
ProdTest( const vector<double> chain )
{
  const Measurement original  = Poisson::Minos( chain.front() );
  const Measurement finalpass = Poisson::Minos( chain.back() );

  vector<Measurement> prodlist;
  prodlist.push_back( original );

  std::stringstream ss;
  ss << chain.front() << flush;

  for( size_t i = 0; i < chain.size()-1; ++i ){
    ss << usr::fstr( "->%lg", chain[i+1] ) << flush;
    prodlist.push_back( Efficiency::Minos( chain[i+1], chain[i] ) );
  }

  const Measurement product = ProdUncorrelated( prodlist );

  cout << fstr( "%30s & %15s & %15s & %15s",
    ss.str(),
    fmt::decimal( product.CentralValue(),  8 ),
    fmt::decimal( product.AbsUpperError(), 8 ),
    fmt::decimal( product.AbsLowerError(), 8 ) )
       << endl;
}

int
main( int argc, char const* argv[] )
{
  cout << separator() << endl
       << ">>>> Addition testing: Sum to 100" << endl;
  {
    Measurement a = Poisson::Minos( 10 );
    Measurement b = Poisson::Minos( 20 );
    Measurement c = Poisson::Minos( 30 );
    Measurement d = Poisson::Minos( 40 );
    Measurement e = Poisson::Minos( 50 );
    Measurement f = Poisson::Minos( 60 );
    Measurement g = Poisson::Minos( 70 );
    Measurement h = Poisson::Minos( 80 );
    Measurement i = Poisson::Minos( 90 );

    cout << "Direct value: " << fmt::decimal( Poisson::Minos( 100 ), 4 ) << endl;
    AddTest( {e, e} );
    AddTest( {d, f} );
    AddTest( {c, g} );
    AddTest( {b, h} );
    AddTest( {a, i} );
    AddTest( {a, c, c, c} );
    AddTest( {b, b, b, b, b} );
    AddTest( {a, a, a, a, a, a, a, a, a, a} );

    cout << endl;
  }

  cout << separator() << endl
       << ">>>> Addition testing: Sum to 1000" << endl;
  {
    Measurement a = Poisson::Minos( 100 );
    Measurement b = Poisson::Minos( 200 );
    Measurement c = Poisson::Minos( 300 );
    Measurement d = Poisson::Minos( 400 );
    Measurement e = Poisson::Minos( 500 );
    Measurement f = Poisson::Minos( 600 );
    Measurement g = Poisson::Minos( 700 );
    Measurement h = Poisson::Minos( 800 );
    Measurement i = Poisson::Minos( 900 );

    cout << "Direct value: " << fmt::decimal( Poisson::Minos( 1000 ), 4 ) << endl;
    AddTest( {e, e} );
    AddTest( {d, f} );
    AddTest( {c, g} );
    AddTest( {b, h} );
    AddTest( {a, i} );
    AddTest( {a, c, c, c} );
    AddTest( {b, b, b, b, b} );
    AddTest( {a, a, a, a, a, a, a, a, a, a} );

    cout << endl;
  }


  cout << separator() << endl
       << ">>>> Addition testing: Partial addition to 100" << endl;
  {
    Measurement a = Poisson::Minos( 20 );
    const std::string test( "%30s & %25s\n" );
    usr::fout( test, "Direct value: ",
      fmt::decimal( Poisson::Minos( 100 ), 4 ) );
    usr::fout( test, "(20+20+20+20+20)",
      fmt::decimal( Sum( a, a, a, a, a ), 4  ) );
    usr::fout( test, "(20+20)+(20+20+20)",
      fmt::decimal( Sum( Sum( a, a ), Sum( a, a, a ) ), 4  ) );
    usr::fout( test, "((20+20)+(20+20)+20)",
      fmt::decimal( Sum( Sum( a, a ), Sum( a, a ), a ), 4  ) );
    usr::fout( test, "(((20+20)+20)+20)+20)",
      fmt::decimal( Sum( Sum( Sum( Sum( a, a ), a ), a ), a ), 4  ) );
  }


  cout << separator() << endl
       << ">>>> Product testing: Product to 50" << endl;
  {
    cout << "Direct value : "
         << fmt::decimal( Poisson::Minos( 50 ), 4 ) << endl;
    ProdTest( {100, 50} );
    ProdTest( {100, 75, 50} );
    ProdTest( {100, 75, 60, 50} );
    ProdTest( {200, 100, 75, 60, 50} );
    cout << endl;
  }

  cout << separator() << endl
       << ">>> Product test: Partial product " << endl;
  {
    Measurement a  = Poisson::Minos( 200 );
    Measurement e1 = Efficiency::Minos( 100, 200 );
    Measurement e2 = Efficiency::Minos( 75, 100 );
    Measurement e3 = Efficiency::Minos( 60, 75 );
    Measurement e4 = Efficiency::Minos( 50, 60 );

    cout << "Direct value: " << fmt::decimal( Poisson::Minos( 50 ), 4 ) << endl;
    cout << "(a*b*c*d*e)" << fmt::decimal( Prod( a, e1, e2, e3, e4 ), 4  ) << endl;
    cout << "(a*b)*(c*d*e)" << fmt::decimal( Prod( Prod( a, e1 ), Prod( e2, e3, e4 ) ), 4  ) << endl;
    cout << "(a*b)*(c*d)*e" << fmt::decimal( Prod( Prod( a, e1 ), Prod( e2, e3 ), e4 ), 4  ) << endl;
    cout << "((((a*b)*c)*d)*e)" << fmt::decimal( Prod( Prod( Prod( Prod( a, e1 ), e2 ), e3 ), e4 ), 4  ) << endl;
  }

  cout << separator() << endl
       << ">>> Another product test: real scale factor tests" << endl;
  {
    Measurement s( 1, 0.00001, 0.00001 );
    Measurement a( 0.9867, 0.0060, 0.0060 );
    Measurement b( 0.9791, 0.0132, 0.0132 );
    Measurement c( 0.9384, 0.0022, 0.0023 );

    cout << fmt::decimal( a, 3 ) << endl;
    cout << fmt::decimal( b, 3 ) << endl;
    cout << fmt::decimal( c, 3 ) << endl;
    cout << fmt::decimal( Prod( a, b, c ), 3 ) << endl;
    cout << fmt::decimal( s *= a, 3 ) << endl;
    cout << fmt::decimal( s *= b, 3 ) << endl;
    cout << fmt::decimal( s *= c, 3 ) << endl;
  }

  cout << separator() << endl
       << ">> Specical case testing " << endl;
  {
    Measurement a( 11394.14521, +3455.26185, +2904.08486 );
    Measurement b( 1.65585, +1.88848, +0.83634 );

    cout << fmt::decimal( a, 5 ) <<endl;
    cout << fmt::decimal( b, 5 ) <<endl;
    a += b;
    cout << fmt::decimal( a, 5 ) << endl;
  }

  cout << separator() << endl
       << "\n>>> Scale case testing" << endl;
  {
    Measurement a( 0, 0, 0 );
    Measurement b = Poisson::Minos( 10 );

    cout << fmt::decimal( b, 5 ) << endl;
    cout << fmt::decimal( a+b, 5 ) << endl;
    cout << fmt::decimal( b+a, 5 ) << endl;
  }

  cout << separator() << endl
       << ">>> Product test" << endl;
  {
    Measurement a( 1, 0.00001, 0.00001 );
    Measurement b( 0.989, 0.0015, 0.0015 );
    Measurement c( 0.989, 0.0015, 0.0015 );
    Measurement d( 0.989, 0.0015, 0.0015 );

    cout << fmt::decimal( Prod( a, b, c ), 5 ) << endl;
    cout << fmt::decimal( Prod( a, b ), 5 ) << endl;
    cout << fmt::decimal( Prod( a, c ), 5 ) << endl;
    cout << fmt::decimal( Prod( b, c ), 5 ) << endl;
    cout << fmt::decimal( Prod( b, c, d ), 5 ) << endl;
  }

  return 0;
}
