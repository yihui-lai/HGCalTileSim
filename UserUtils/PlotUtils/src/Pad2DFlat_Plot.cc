#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/Pad2DFlat.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/Pad2DFlat.hpp"
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "CmdSetAttr.hpp"
#include "TLegendEntry.h"

namespace usr {

namespace plt {

/**
 * Plotting a two dimensional histogram support the following options:
 *
 * - Plot2DF: Defining how the the data should be represented. The following
 *   types are supported.
 *   - `plot2df::heat`: Heat map, corresponds to the `"COLZ"` options for TH2
 *     plotting in @ROOT
 *   - `plot2df::cont` Plotting the contour lines, using the "CONT3" options for
 *      TH2 plotting
 *   - `plot2df::heatcont`: Plotting both the heat map and contour lines.
 *
 * - EntryText: Adding text for display in the legend. Currently only the `cont`
 *   and `heatcont` will add a legend entry with the line styling displayed. The
 *   `heat` option assumes that this object will be used for the z axis colour
 *   pallet, so it will only add a text entry in the legend without a display.
 */
TH2D&
Pad2DFlat::PlotHist( TH2D& hist, const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args( arglist,
                              { Plot2DF( plot2df::heat )  } );

  const int opt = args.Get<Plot2DF>();

  const std::string legopt = opt == plot2df::heat ? "" :
                             opt == plot2df::heatcont ? "L" :
                             opt == plot2df::cont ? "L" :
                             "PFLE";

  for( const auto&& func : *( hist.GetListOfFunctions() ) ){
    func->SetBit( TF1::kNotDraw, true );
  }

  // Forcing stuff to not be displayed
  hist.SetTitle( "" );
  hist.SetStats( 0 );

  // Parsing plotting flag
  switch( opt ){
  case  plot2df::heat:
    PlotObj( hist, "COLZ" );
    SetAxisFont();
    break;
  case plot2df::heattext:
    PlotObj( hist, "COLZ TEXT" );
    SetAxisFont();
    break;
  case   plot2df::heatcont:
    PlotObj( hist, "COLZ" );
    SetAxisFont();
    PlotObj( hist, "CONT3 SAME" );
    break;
  case  plot2df::cont:
    PlotObj( hist, "CONT3 SAME" );
    SetAxisFont();
    break;
  case plot2df_dummy:
    PlotObj( hist, ( args.Get<Plot2DF>().str() + " SAME" ).c_str() );
    break;

  default:
    std::cerr << "Skipping over invalid value" << std::endl;
    break;
  }

  if( args.Has<EntryText>() ){
    _legend.AddEntry( &hist, args.Get<EntryText>().c_str(), legopt.c_str() );
  }

  SetLineAttr( hist, args );
  SetFillAttr( hist, args );
  SetMarkAttr( hist, args );

  return hist;
}

/**
 * Plotting a two dimensional graph support the following options:
 *
 * - Plot2DF: Defining how the the data should be represented. The following
 *   types are supported.
 *   - `plot2df::heat`: Heat map, corresponds to the "CONT4" options for
 *     TGraph2D plotting in @ROOT
 *   - `plot2df::cont` Plotting the contour lines, using the "CONT3" options for
 *      TH2 plotting
 *   - `plot2df::heatcont`: Plotting both the heat map and contour lines.
 *
 * - EntryText: Adding text for display in the legend. Currently only the `cont`
 *   and `heatcont` will add a legend entry with the line styling displayed. The
 *   `heat` option assumes that this object will be used for the z axis colour
 *   pallet, so it will only add a text entry in the legend without a display.
 * @brief
 */
TGraph2D&
Pad2DFlat::PlotGraph( TGraph2D& graph, const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args( arglist,
                              { Plot2DF( plot2df::heat )  } );

  const int opt = args.Get<Plot2DF>();

  const std::string legopt = opt == plot2df::heat ? "" :
                             opt == plot2df::heatcont ? "L" :
                             opt == plot2df::cont ? "L" :
                             "PFLE";

  if( !GetAxisObject() ){
    auto& axishist = MakeObj<TH2D>(
      ( "axishist" + RandomString( 12 ) ).c_str(), "",
      10, graph.GetXminE(), graph.GetXmaxE(),
      10, graph.GetYminE(), graph.GetYmaxE()
      );
    axishist.SetStats( 0 );
    PadBase::PlotObj( axishist, "AXIS" );
    SetAxisFont();
  }

  for( const auto&& func : *( graph.GetListOfFunctions() ) ){
    func->SetBit( TF1::kNotDraw, true );
  }

  switch( opt ){
  case plot2df::heat:
    PlotObj( graph, "CONT4Z SAME" );
    break;
  case  plot2df::heatcont:
    PlotObj( graph, "CONT4Z SAME" );
    PlotObj( graph, "CONT3 SAME" );
    break;
  case plot2df::cont:
    PlotObj( graph, "CONT3 SAME" );
    break;
  case plot2df_dummy:
    PlotObj( graph, args.Get<Plot2DF>().c_str() );
    break;

  default:
    std::cerr << "Skipping over invalid value" << std::endl;
    break;
  }

  if( args.Has<EntryText>() ){
    _legend.AddEntry( &graph, args.Get<EntryText>().c_str(), legopt.c_str() );
  }

  SetLineAttr( graph, args );
  SetFillAttr( graph, args );
  SetMarkAttr( graph, args );

  return graph;
}

/**
 * Plotting of a 2 dimensional graph is done by generating a TGraph2D object via
 * the sampling the function. All other plotting functionalities are identical
 * to the the PlotGraph function.
 */
TGraph2D&
Pad2DFlat::PlotFunc( TF2& func, const std::vector<RooCmdArg>& arglist )
{
  static const unsigned sep = 300;

  const double xmin = func.GetXmin();
  const double xmax = func.GetXmax();
  const double ymin = func.GetYmin();
  const double ymax = func.GetYmax();
  const double dx   = ( xmax-xmin )/sep;
  const double dy   = ( ymax-ymin )/sep;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  for( unsigned i = 0; i < sep; ++i ){
    for( unsigned j = 0; j < sep; ++j ){
      x.push_back( xmin + ( i + 0.5 ) * dx );
      y.push_back( ymin + ( j + 0.5 ) * dy );
      z.push_back( func.Eval( x.back(), y.back() ) );
    }
  }

  TGraph2D& graph = MakeObj<TGraph2D>( x.size(),
    x.data(), y.data(), z.data() );

  graph.SetName( usr::RandomString( 12 ).c_str() );

  return PlotGraph( graph, arglist );
}

/**
 * Plotting a regular 1D graph object (TGraph as oppose to TGraph2D) object
 * supports the following commands:
 * - PlotType: specifying how the data should be represented. The supported plot
 *   types are:
 *   - `plottype::scatter`: Scatter point with error bars (no end ticks)
 *   - `plottype::simplefunc`: Simple polyline joining the points.
 * - EntryText: Adding text into the legend.
 */
TGraph&
Pad2DFlat::Plot1DGraph( TGraph& graph, const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args( arglist,
    {PlotType(plottype::scatter) } );

  const int opt =args.Get<PlotType>();

  std::string legopt = "PLE";

  // Parsing plotting flag
  if( opt == plottype::scatter ){
    PadBase::PlotObj( graph, "PZ0" );
  } else if( opt == plottype::simplefunc ){
    PadBase::PlotObj( graph, "L" );
    legopt = "L";
  } else if( opt == plottype_dummy ){
    // Special case for raw options parsing. (Must remove axis and add string)
    PadBase::PlotObj( graph, args.Get<PlotType>().c_str() );
  } else {// Skipping over stuff
    std::cerr << "Skipping over invalid value" << std::endl;
  }

  // Resetting underlying histogram range:
  graph.GetHistogram()->SetMaximum( usr::plt::GetYmax( graph ) );
  graph.GetHistogram()->SetMinimum( usr::plt::GetYmin( graph ) );

  if( args.Has<EntryText>() ){
    _legend.AddEntry( &graph, args.Get<EntryText>().c_str(), legopt.c_str() );
  }

  SetLineAttr( graph, args );
  SetFillAttr( graph, args );
  SetMarkAttr( graph, args );

  return graph;
}

/**
 * In 2D plotting, we are assuming the important stuff is being plotted first:
 * Important takes up more area (with heat plots, for example), so we will not
 * need to reverse the legend entry as with the Pad1D objects. The TLegend will
 * be generated on the left side of the frame. The size will be adjusted to have
 * everything fit snugly in the legend.
 */
void
Pad2DFlat::MakeLegend()
{
  // Early exit if Legend wasn't requested
  if( !_legend.GetListOfPrimitives() ){ return; }
  if( !_legend.GetListOfPrimitives()->GetEntries() ){ return; }
  TPad::cd();
  float width  = 0;
  float height = 1.2*LineHeight() * _legend.GetListOfPrimitives()->GetEntries();

  for( const auto&& obj : *_legend.GetListOfPrimitives() ){
    const char* label = ( (TLegendEntry*)obj )->GetLabel();
    TLatex* textmp    = new TLatex( 0, 0, label );
    textmp->SetTextFont( FontFace() );
    textmp->SetTextSize( FontSize() );
    width = std::max( int(width),
      int(1.3 *textmp->GetXsize() * AbsWidth() + 5) );
    delete textmp;
  }

  width *= 1.1;// Relieving spacing a little
  width += 1.0 * LineHeight();// Reserving space for legend icon boxes.

  const float xmin = 0.5 * float( LineHeight() ) / AbsWidth();
  const float ymax = 1   - GetTopMargin();
  const float xmax = GetLeftMargin() - 5*float( LineHeight() ) / AbsWidth();
  const float ymin = ymax - height/ AbsHeight();
  _legend.SetX1NDC( xmin );
  _legend.SetX2NDC( xmax );
  _legend.SetY1NDC( ymin );
  _legend.SetY2NDC( ymax );

  if( !TPad::FindObject( &_legend ) ){
    PadBase::PlotObj( _legend );
  }
}

}/** plt */

}/* usr */
