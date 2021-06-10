/**
 * @file Pad2DFlat.hpp
 * @author your name (you@domain.com)
 * @brief
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/PlotUtils/interface/Pad2DFlat.hpp"
#else
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/PlotUtils/Pad2DFlat.hpp"
#endif


namespace usr {

namespace plt {

Pad2DFlat::Pad2DFlat( const PadSize& padsize ):
  PadBase(padsize)
{
  init_legend();
}

Pad2DFlat::Pad2DFlat( const PadSize& padsize ,
  const double xmin, const double xmax,
  const double ymin, const double ymax ):
  PadBase(padsize),
  _frame(xmin,xmax,ymin,ymax)
{
  init_legend();
}


Pad2DFlat::~Pad2DFlat()
{}

void Pad2DFlat::init_legend()
{
  _legend.SetX1( 0 ) ;
  _legend.SetX2( GetLeftMargin() );
  _legend.SetY1( GetBottomMargin() );
  _legend.SetY2( 1- GetTopMargin() );
  _legend.SetTextFont( ParentCanvas().Font().fontface() );
  _legend.SetTextSize( ParentCanvas().Font().size() );
  _legend.SetBorderSize( 0 );
  _legend.SetFillColorAlpha( 0, 0 );
}

void Pad2DFlat::InitDraw()
{
  // Early exit if no RooPlot is required;
  if( !_frame.AxisHistPtr() ){ return; }

  // Drawing the frame histogram for defining axis
  PadBase::PlotObj( _frame.AxisHist(), "AXIS" );
  // Setting up the name and title of RooPlot:
  SetAxisFont();
}

void Pad2DFlat::Finalize()
{
  TPad::cd();
  PadBase::Finalize();
  MakeLegend();
  MakeLegend();  // Must run twice... not sure why
}

}

}
