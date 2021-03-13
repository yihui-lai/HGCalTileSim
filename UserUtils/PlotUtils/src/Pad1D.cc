/**
 * @file    Pad1D.cc
 * @brief   Implementation of the Pad1D member functions
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/PlotUtils/interface/Pad1D.hpp"
#else
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/PlotUtils/Pad1D.hpp"
#endif

#include "TGaxis.h"
#include <limits>

namespace usr  {

namespace plt {

/**
 * @brief Constructing a Pad1D using a @ROOT{RooRealVar}
 *
 * @details This will make the Pad1D instance properly spawn the _frame that is
 * capable for generating plotting objects for RooFit data objects.
 */
Pad1D::Pad1D( const PadSize& size, const RangeByVar& range ) :
  PadBase( size ),
  _frame( range ),
  _workingstack( nullptr ),
  _datamax( 1 ),
  _datamin( 0.3 ),
  _prevrangetype( rangetype::aut ),
  _legendposition( align::top_right )
{
  // Common setup for pads
  TPad::SetTicks( 1, 1 );
  _init_legend();
}

/**
 * @brief Construct a Pad1D with a limit x range.
 */
Pad1D::Pad1D( const PadSize& size, const double min, const double max ) :
  PadBase( size ),
  _frame( min, max ),
  _workingstack( nullptr ),
  _datamax( 1 ),
  _datamin( 0.3 ),
  _prevrangetype( rangetype::aut ),
  _legendposition( align::top_right )
{
  // Common setup for pads
  TPad::SetTicks( 1, 1 );
  _init_legend();
}

/**
 * @brief Construct a Pad1D with only dimension specifications
 *
 * This would the be most typical construction for the plotting of ROOT objects.
 * Note that the axis range would be determined by the object that is plotted
 * first into the Pad. So take care of plotting order if that is an issue..
 */
Pad1D::Pad1D( const PadSize& size ) :
  PadBase( size ),
  _workingstack( nullptr ),
  _datamax( 1 ),
  _datamin( 0.3 ),
  _prevrangetype( rangetype::aut ),
  _legendposition( align::top_right )
{
  // Common setup for pads
  TPad::SetTicks( 1, 1 );
  _init_legend();
}

/**
 * @brief Additional routines after pad spawning
 *
 * In case a RooPlot is requested, ust that to draw the frame first to
 * avoid graphs messing up the plot range.
 */
void
Pad1D::InitDraw()
{
  // Setting axis digits
  TGaxis::SetMaxDigits( 4 );

  // Early exit if no RooPlot is required;
  if( !_frame.AxisHistPtr() ){ return; }

  // Drawing the frame histogram for defining axis
  PadBase::PlotObj( _frame.AxisHist(), "AXIS" );
  // Setting up the name and title of RooPlot:
  SetAxisFont();

  if( _frame.getPlotVar() ){
    SetHistAxisTitles(
      _frame.getPlotVar()->GetTitle(),
      _frame.getPlotVar()->getUnit()
      );
  }
}

/**
 * @brief additional steps before canvas saving.
 *
 * Additional steps:
 * 1. Redrawing axis to ensure nothing is covering the frame.
 * 2. Making the legend object, and ensuring all specificed legend entries
 *    are placed in the legend.
 */
void
Pad1D::Finalize()
{
  TPad::cd();
  PadBase::Finalize();
  FinalizeLegend( _legendposition );
  FixVLines();
  TPad::RedrawAxis();
}

/**
 * Since all ownership is handled by _frame, do additional deletion class is
 * needed.
 */
Pad1D::~Pad1D()
{}

/**
 * Aside from calling the TPad::SetLogy() function, the AutoSetYRange is also
 * called to adjust the Y axis range appropriately.
 */
void
Pad1D::SetLogy( int flag )
{
  TPad::SetLogy( flag );
  if( flag && !CheckLogy() ){
    std::cerr << "[PAD1D] Data-objects y-values in Pad1D is not "
      "positive definite. Plot may break!" << std::endl;
    // TPad::SetLogy( 0 );
  }
  AutoSetYRange();
}

/**
 * @brief Simple passthrough function
 */
void
Pad1D::SetLogx( int flag )
{
  TPad::SetLogx( flag );
  if( flag && !CheckLogx() ){
    std::cerr << "[PAD1D] Data-objects x-values in Pad1D is not "
      "positive definite. Unsetting the log option." << std::endl;
    TPad::SetLogx( 0 );
  }
}

/**
 * @brief Checking if the data in the plot objects are positive definite.
 *
 * To avoid the case where there are negative numbers or zeros in the plot data
 * when setting the plot range to be log scaled. All plot objects with values
 * less than 0 with have it's contents shifted to the minimal positive double.
 */
bool
Pad1D::CheckLogy() const
{
  return GetDataMin() > 0 ;
}

bool
Pad1D::CheckLogx() const
{
  return Xaxis().GetXmin() >= 0;
}

}/* plt */

}/* usr  */
