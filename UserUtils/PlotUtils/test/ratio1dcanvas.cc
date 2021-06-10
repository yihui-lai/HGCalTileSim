/**
 * @file    ratiocanvas.cc
 * @brief   Testing the default value settings of the ratiocanvas class
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#include "UserUtils/MathUtils/interface/RooFitExt.hpp"
#else
#include "UserUtils/PlotUtils/Ratio1DCanvas.hpp"
#include "UserUtils/MathUtils/RooFitExt.hpp"
#endif

#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TH1D.h"
#include "TRandom3.h"
#include <iostream>

namespace plt = usr::plt;

template<typename T>
void
SetDataStyle( T& x )
{
  x.SetMarkerStyle( 20 );
  x.SetMarkerSize(    0.5 );
  x.SetMarkerSize( kBlack );
}

int
main( int argc, char* argv[] )
{
  {// Root object testing.
    TH1D hist1( "h1", "", 21, -10.5, 10.5 );
    TH1D hist2( "h2", "", 21, -10.5, 10.5 );
    TH1D hist3( "h3", "", 21, -10.5, 10.5 );
    TH1D histsum( "sig", "", 21, -10.5, 10.5 );
    TH1D data( "dat", "", 21, -10.5, 10.5 );
    TRandom3 rnd;

    for( size_t i = 0; i < 5000; i++ ){
      hist1.Fill( rnd.Gaus( -2, 2 ), 0.3 );
      hist2.Fill( rnd.Gaus( -2.3, 2.5 ), 0.3*rnd.Uniform( 0.0, 0.6 ) );
      hist3.Fill( rnd.Gaus( -1.7, 3.0 ), 0.3*rnd.Uniform( 0.0, 0.4 ) );
      data.Fill( rnd.Gaus( -2, 2.2 ), 0.5 );
    }

    histsum.Add( &hist1 );
    histsum.Add( &hist2 );
    histsum.Add( &hist3 );

    // Artificially magnifying uncertainty
    for( int i = 0; i < histsum.GetNcells(); i++ ){
      histsum.SetBinError( i,
        std::min( histsum.GetBinError( i )*1.5*rnd.Uniform( 0.8, 1.2 ),
          histsum.GetBinContent( i ) ) );
      data.SetBinError( i,
        std::min( histsum.GetBinError( i )*2*rnd.Uniform( 0.8, 1.2 ),
          histsum.GetBinContent( i ) ) );
    }

    plt::Ratio1DCanvas c;

    c.PlotHist( hist3,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg_{1}: VV" ),
      plt::FillColor( plt::col::blue ),
      plt::LineColor( plt::col::blue ) );
    c.PlotHist( hist2,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg_{2}: HH" ),
      plt::FillColor( plt::col::red ),
      plt::LineColor( plt::col::red ) );
    c.PlotHist( hist1,
      plt::PlotType( plt::histstack ),
      plt::EntryText( "Bkg_{3}: TT" ),
      plt::FillColor( plt::col::cyan ),
      plt::LineColor( plt::col::cyan ) );
    c.PlotHist( histsum,
      plt::PlotType( plt::histerr ),
      plt::EntryText( "Bkg_{1} Unc." ),
      plt::FillStyle( plt::sty::filldotdense ),
      plt::FillColor( plt::col::darkgray ) );
    c.PlotHist( data,
      plt::PlotType( plt::scatter ),
      plt::EntryText( "Data" ),
      plt::MarkerStyle( 20 ),
      plt::MarkerSize( 0.2 ),
      plt::LineColor( plt::col::black ) );

    c.PlotScale( histsum, histsum,
      plt::PlotType( plt::histerr ),
      plt::TrackY( plt::TrackY::none ) );
    c.PlotScale( data, histsum,
      plt::PlotType( plt::scatter ) );

    c.TopPad().SetHistAxisTitles( "M_{THTHTH}", plt::unit::GeVcc );
    c.BottomPad().SetHistAxisTitles( "M_{THTHTH}", plt::unit::GeVcc );
    c.BottomPad().Yaxis().SetTitle( "Data/Sim." );

    c.DrawCMSLabel( "Ratio1DCanvas", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SaveAsPNG(         "image/ratio1dcanvas.png",  72 );
    c.SaveAsPNG( "image/ratio1dcanvas_highres.png", 300 );
    c.SaveAsPDF( "image/ratio1dcanvas_highres.pdf" );
  }
  {// RooFit Testing
    RooRealVar x( "x", "x", -10, 10, plt::unit::GeVcc.c_str() );
    RooRealVar m( "m", "m", -2, -5, 5 );
    RooRealVar s( "s", "s", 1.5, 0, 5 );
    RooGaussian g( "g", "g", x, m, s );
    RooGaussian gx( "gx", "gx", x,
                    RooFit::RooConst( -1.85 ), RooFit::RooConst( 1.40 ) );
    RooDataSet* d     = g.generate( RooArgSet( x ), 1000 );

    RooFitResult* fit = usr::FitPDFToData( g, *d, RooFit::Save() );

    plt::Ratio1DCanvas c( x );
    auto& datgraph = c.PlotData( d,
      plt::EntryText( "Fake Data" ),
      plt::MarkerStyle( 20 ),
      plt::MarkerSize( 0.2 ),
      plt::LineColor( kBlack ) );

    auto& altgraph = c.PlotPdf( gx,
      plt::EntryText( "Alt. Bkg. Model" ),
      plt::PlotUnder( datgraph ),
      plt::LineColor( plt::col::red ),
      plt::FillColor( plt::col::pink ) );

    auto& fitgraph = c.PlotPdf( g,
      plt::VisualizeError( *fit, 1, false ),
      plt::EntryText( "Fitted Bkg." ),
      plt::PlotUnder( datgraph ),
      plt::LineColor( usr::plt::col::blue ),
      plt::FillColor( usr::plt::col::cyan )  );

    c.PlotScale( fitgraph, fitgraph,
      plt::TrackY( plt::TrackY::none ), plt::PlotType( plt::fittedfunc ) );
    c.PlotScale( altgraph, fitgraph,
      plt::TrackY( plt::TrackY::none ), plt::PlotType( plt::simplefunc ) );
    c.PlotScale( datgraph, fitgraph,
      plt::PlotType( plt::scatter ) );

    c.BottomPad().Yaxis().SetTitle( "Data/Fit" );

    c.DrawCMSLabel( "Ratio1DCanvas", "CWS" );
    c.DrawLuminosity( 133.7 );
    c.SetLogy(1);
    c.SaveAsPNG( "image/ratio1dcanvas_roofit.png",  72 );
    c.SaveAsPNG( "image/ratio1dcanvas_roofit_highres.png", 300 );
    c.SaveAsPDF( "image/ratio1dcanvas_roofit.pdf" );
  }

  return 0;
}
