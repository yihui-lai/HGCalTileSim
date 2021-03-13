/**
 * @file    flat2canvas.cc
 * @brief   Testing the default value settings of the ratiocanvas class
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#else
#include "UserUtils/PlotUtils/Flat2DCanvas.hpp"
#endif

#include "TF2.h"
#include "TGraph2D.h"
#include "TH2D.h"

int
main( int argc, char* argv[] )
{
  TH2D hist( "hist", "hist", 20, -5, 5, 20, -5, 5 );
  TF2 xyg( "xyg", "xygaus", -5, 5, -5, 5 );
  xyg.SetParameters( 1, 1.5, 2, 1.2, 2 );
  TF2 xyg2( "xyg2", "xygaus", -5, 5, -5, 5 );
  xyg2.SetParameters( 1, 1.2, 1.8, 1.3, 2.2 );
  hist.FillRandom( "xyg", 100000 );

  hist.SetContour( 50 );

  usr::plt::Flat2DCanvas c;

  auto& gh = c.PlotHist( hist,
    usr::plt::Plot2DF( usr::plt::heat ) );
  auto& gf1 = c.PlotFunc( xyg,
    usr::plt::Plot2DF( usr::plt::cont ),
    usr::plt::EntryText( "Model 1" )
    );
  auto& gf2 = c.PlotFunc( xyg2,
    usr::plt::Plot2DF( usr::plt::cont ),
    usr::plt::EntryText( "Model 2" )
    );

  c.DrawCMSLabel("2D plotting", "CWS" );
  c.DrawLuminosity( 133.7 );

  c.Xaxis().SetTitle( "My X value" );
  c.Yaxis().SetTitle( "My y value" );
  c.Zaxis().SetTitle( "Readout" );

  gf1.SetLineColor( kRed );
  gf2.SetLineColor( kBlue );

  c.SaveAsPNG( "image/flat2dcanvas.png" , 72 );
  c.SaveAsPNG( "image/flat2dcanvas_highres.png" , 300 );
  c.SaveAsPDF( "image/flat2dcanvas.pdf" );
}
