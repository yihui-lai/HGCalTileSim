/**
 * @file    simple1dcanvas.cc
 * @brief   Testing the simple 1d Canvas.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#else
#include "UserUtils/PlotUtils/Simple1DCanvas.hpp"
#include "UserUtils/MathUtils/RooFitExt.hpp"
#endif

#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TRandom3.h"

namespace plt = usr::plt;

int
main( int argc, char* argv[] )
{
  {// Histogram testing
    TRandom3 r;
    TH1D h1( "h1", "", 21, -10, 10 );
    TH1D h2( "h2", "", 21, -10, 10 );
    TH1D h3( "h3", "", 21, -10, 10 );

    for( int i = 0; i < 5000; ++i ){
      h1.Fill( r.Gaus( 0, 3 ) );
      h2.Fill( r.Gaus( 1, 4 ) );
      h3.Fill( r.Gaus( -2, 1 ), r.Gaus( 0.5, 0.1 ) );
    }

    plt::Simple1DCanvas c;
    c.PlotHist( h1,
      plt::PlotType( plt::scatter ),
      plt::EntryText( "Bkg. 1" ),
      plt::LineColor( usr::plt::col::black ),
      plt::MarkerStyle( usr::plt::sty::mkrcircle ),
      plt::MarkerColor( usr::plt::col::black ),
      plt::MarkerSize( 0.2 ) );
    c.PlotHist( h2,
      plt::PlotType( plt::hist ),
      plt::EntryText( "Bkg. 2" ),
      plt::LineColor( usr::plt::col::blue ) );
    c.PlotHist( h3,
      plt::PlotType( plt::histerr ),
      plt::EntryText( "Signal" ),
      plt::FillColor( usr::plt::col::darkgray ),
      plt::FillStyle( plt::sty::filldotdense ) );
    c.DrawCMSLabel( "ROOTOBJ TEST", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SetHistAxisTitles( "P_{t}", plt::unit::GeVc );

    c.SetLogy( 1 );
    c.SaveAsPNG( "image/simple1dcanvas_log.png",         72 );
    c.SaveAsPNG( "image/simple1dcanvas_log_highres.png", 300 );
    c.SaveAsPDF( "image/simple1dcanvas_log.pdf" );
  }

  {// Stacking test
    TRandom3 r;
    TH1D h1( "h1_new", "", 21, -10, 10 );
    TH1D h2( "h2_new", "", 21, -10, 10 );
    TH1D h3( "h3_new", "", 21, -10, 10 );

    for( int i = 0; i < 5000; ++i ){
      h1.Fill( r.Gaus( 0, 3 ) );
      h2.Fill( r.Gaus( 1, 4 ) );
      h3.Fill( r.Gaus( -2, 1 ), r.Gaus( 0.5, 0.1 ) );
    }

    plt::Simple1DCanvas c;
    c.PlotHist( h1,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg. 1" ),
      plt::FillStyle( plt::sty::fillsolid ),
      plt::FillColor( plt::col::cyan ),
      plt::LineColor( plt::col::cyan ) );
    c.PlotHist( h2,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg. 2" ),
      plt::FillStyle( plt::sty::filldotdense ),
      plt::FillColor( plt::col::pink ),
      plt::LineColor( plt::col::pink ) );
    c.PlotHist( h3,
      plt::PlotType( plt::hist ),
      plt::EntryText( "Sig" ),
      plt::LineColor( usr::plt::col::blue ) );

    c.DrawCMSLabel( "STACK TEST", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SetHistAxisTitles( "P_{t}", plt::unit::GeVc );

    c.SaveAsPNG( "image/simple1dcanvas_stack.png",         72 );
    c.SaveAsPNG( "image/simple1dcanvas_stack_highres.png", 300 );
    c.SaveAsPDF( "image/simple1dcanvas_stack.pdf" );
    c.SetLogy( 1 );
    c.SaveAsPNG( "image/simple1dcanvas_stack_log.png",         72 );
    c.SaveAsPNG( "image/simple1dcanvas_stack_log_highres.png", 300 );
    c.SaveAsPDF( "image/simple1dcanvas_stack_log.pdf" );
  }

  // ROOFIT testing
  {
    RooRealVar x( "x", "x", -10, 10, "GeV/#it{c}^2" );
    RooRealVar m( "m", "m", 0.3, -5, 5 );
    RooRealVar s( "s", "s", 2, 0, 5 );
    RooGaussian g( "g", "g", x, m, s );
    RooDataSet* d     = g.generate( RooArgSet( x ), 800 );
    RooFitResult* fit = usr::FitPDFToData( g, *d, RooFit::Save() );

    plt::Simple1DCanvas c2( x );

    auto& dgraph = c2.PlotData( d,
      plt::EntryText("Data?"),
      plt::MarkerStyle( plt::sty::mkrcircle ),
      plt::MarkerSize( 0.5 ),
      plt::MarkerColor( plt::col::black ),
      plt::LineColor( plt::col::black ) );

    c2.PlotPdf( g,
      RooFit::VisualizeError( *fit, 1, false ),
      plt::EntryText( "My Fit" ),
      plt::PlotUnder( dgraph ),
      plt::LineColor( plt::col::blue ),
      plt::FillColor( plt::col::cyan )
      );

    c2.DrawCMSLabel( "ROOFIT test", "CWS" );
    c2.DrawLuminosity( 133.7 );
    c2.SetHistAxisTitles( "M_{HHHHH}", plt::unit::GeVcc );

    c2.SaveAsPNG( "image/simple1dcanvas_roofit.png",         72 );
    c2.SaveAsPNG( "image/simple1dcanvas_roofit_highres.png", 300 );
    c2.SaveAsPDF( "image/simple1dcanvas_roofit.pdf" );
  }

  return 0;
}
