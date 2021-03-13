/**
 * @file    Pad1D_Axis.cc
 * @brief   Implementation of axis access and range adjustment function
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/PlotUtils/interface/Pad1D.hpp"
#else
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/PlotUtils/Pad1D.hpp"
#endif

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

namespace usr  {

namespace plt  {

/**
 * A static constant for getting defining something as not-found
 */
static double unfound = std::numeric_limits<double>::quiet_NaN();

/**
 * @brief Returning the TObject used for drawing the axis frame.
 *
 * Getting the pointer to the TObject responsible for drawing the axis.
 * Based on how the Various Draw functions works, this will be the first object
 * on the TPad that is either a TH1 or a TGraph. The code for iterating through
 * the objects in the pad is copied from the  [TPad::RedrawAxis()](https://root.cern.ch/doc/master/TPad_8cxx_source.html#l05243) code.
 */
TObject*
Pad1D::GetAxisObject() const
{
  if( !GetListOfPrimitives() ){ return nullptr; }

  for( const auto&& obj : *GetListOfPrimitives() ){
    if( obj->InheritsFrom( TH1::Class() ) ){
      return obj;
    }
  }

  return nullptr;
}

/**
 * @brief Setting the axis fonts according to canvas settings.
 *
 * Setting the axis font based on the Canvas settings. The axis title should be
 * the same size as the Canvas font, while the label (numbers) would be the
 * scale to the latex "small" format.
 */
void
Pad1D::SetAxisFont()
{
  Xaxis().SetLabelFont( FontFace() );
  Xaxis().SetTitleFont( FontFace() );
  Xaxis().SetTitleSize( FontSize() );
  Xaxis().SetLabelSize( Font().small() );

  Yaxis().SetLabelFont( FontFace() );
  Yaxis().SetTitleFont( FontFace() );
  Yaxis().SetTitleSize( FontSize() );
  Yaxis().SetLabelSize( Font().small() );
}


#define AXISOBJ_ACCESS( ACTION, DEFAULT )                 \
  TObject* axisobj = GetAxisObject();                     \
  if( !axisobj ){                                         \
    return DEFAULT;                                       \
  } else if( axisobj->InheritsFrom( TH1::Class() ) ){     \
    return dynamic_cast<TH1*>( axisobj )->ACTION;         \
  } else if( axisobj->InheritsFrom( TGraph::Class() ) ){  \
    return dynamic_cast<TGraph*>( axisobj )->ACTION;      \
  } else if( axisobj->InheritsFrom( THStack::Class() ) ){ \
    return dynamic_cast<THStack*>( axisobj )->ACTION;     \
  } else {                                                \
    return DEFAULT;                                       \
  }

#define AXISOBJ_ADJUST( ACTION )                          \
  TObject* axisobj = GetAxisObject();                     \
  if( !axisobj ){                                         \
    return;                                               \
  } else if( axisobj->InheritsFrom( TH1::Class() ) ){     \
    return dynamic_cast<TH1*>( axisobj )->ACTION;         \
  } else if( axisobj->InheritsFrom( TGraph::Class() ) ){  \
    return dynamic_cast<TGraph*>( axisobj )->ACTION;      \
  } else if( axisobj->InheritsFrom( THStack::Class() ) ){ \
    return dynamic_cast<THStack*>( axisobj )->ACTION;     \
  }


/**
 * @{
 * @brief returning pointer to axis object.
 * @details
 * ROOT style pointer interface for getting the axis object used for plotting
 * the axis frame.
 */
TAxis*
Pad1D::GetXaxis() const { AXISOBJ_ACCESS( GetXaxis(), nullptr ); }

TAxis*
Pad1D::GetYaxis() const { AXISOBJ_ACCESS( GetYaxis(), nullptr ); }
/** @} */

/**
 * @{
 * @brief   reference interface for getting axis (No null pointer handling)
 */
TAxis&
Pad1D::Xaxis(){ return *GetXaxis(); }

TAxis&
Pad1D::Yaxis(){ return *GetYaxis(); }

const TAxis&
Pad1D::Xaxis() const { return *GetXaxis(); }

const TAxis&
Pad1D::Yaxis() const { return *GetYaxis(); }
/** @} */


/**
 * @{
 * @brief Getting the x axis limits.
 *
 * Since this would not be altered thoughout the plotting function calls, we
 * directly use the in-built function of the TAxis object.
 */
double
Pad1D::GetXaxisMin() const { return Xaxis().GetXmin(); }

double
Pad1D::GetXaxisMax() const { return Xaxis().GetXmax(); }
/** @} */

/**
 * @{
 * @brief Getting y axis limit.
 *
 * The Y axis plotting range is not defined by the Y axis object itself, but
 * rather the GetMaximum() function of the object used for drawing the Y axis
 * of the canvas. This function would return an NAN the the GetAxisObject is
 * not found.
 */
double
Pad1D::GetYaxisMax() const { AXISOBJ_ACCESS( GetMaximum(), unfound ); }

double
Pad1D::GetYaxisMin() const { AXISOBJ_ACCESS( GetMinimum(), unfound ); }
/** @} */

/**
 * @{
 * @brief Setting y axis limit
 *
 * Again, setting the y axis limit should call the SetMaximum() function of the
 * object used for drawing the Y axis.
 */
void
Pad1D::SetYaxisMax( const double x ){ AXISOBJ_ADJUST( SetMaximum( x ) ); }

void
Pad1D::SetYaxisMin( const double x ){ AXISOBJ_ADJUST( SetMinimum( x ) ); }
/** @} */


/**
 * @{
 * @brief getting tracked y limit of data.
 *
 * Directly returning the double stored by the class.
 */
double
Pad1D::GetDataMax() const { return _datamax; }

double
Pad1D::GetDataMin() const { return _datamin; }
/** @} */

/**
 * @{
 * @brief Manually overriding the stored y limit of data
 *
 * Directly returning the double stored by the class.
 */
void
Pad1D::SetDataMax( const double x ){ _datamax = x; }

void
Pad1D::SetDataMin( const double x ){ _datamin = x; }
/** @} */


/**
 * @brief Automatically adjusting y axis range according to tracked min/max
 *        values.
 *
 * Automatically adjusting the Y axis range according the the data's min max
 * value that are tracked by the Pad1D. A few types are present:
 * - rangetype::hist: The Y range optimized for the presentation histograms.
 *   Mainly adjusting the max value, leaving a larger margin above the maximum
 *   peak. See AutoSetYRangeHist().
 * - rangetype::graph: The Y range is optimized for the presentation of graphs.
 *   Adjusting both the min and max range so that the data would fit
 *   comfortably in the frame. See AutoSetYRangeGraph().
 * - rangetype::ratio: The Y range would be symmetric around 1. See
 *   AutoSetYRangeRatio().
 * - rangetype::pull: The Y range would be symmetric around 0. See
 *   AutoSetYRangePull()
 * More details would be given in the implementations in the protected member
 * functions.
 *
 * In the case the the user want's this to be automatically determined, if
 * any histogram like data exists in the pad, this method used would be
 * the rangetype::hist, otherwise if only TGraph-like data exists in the pad,
 * the rangetype::graph would be used.
 */
void
Pad1D::AutoSetYRange( const rangetype type )
{
  rangetype optype =
    type != rangetype::aut ? type :
    _prevrangetype != rangetype::aut ? _prevrangetype :
    rangetype::aut;

  if( optype == rangetype::aut ){
    if( !GetListOfPrimitives() ){
      optype = rangetype::aut;
    }

    for( const auto&& obj : *GetListOfPrimitives() ){
      if( ( obj->InheritsFrom( TH1::Class() )
            && ( (TH1*)obj )->GetXaxis() != _frame.GetXaxis() )
          || ( obj->InheritsFrom( THStack::Class() ) )
           ){
        optype = rangetype::hist;
        break;
      } else if( obj->InheritsFrom( TGraph::Class() ) ){
        optype = rangetype::graph;
      }
    }
  }

  if( optype == rangetype::hist ){
    AutoSetYRangeHist();
  } else if( optype == rangetype::graph ){
    AutoSetYRangeGraph();
  } else if( optype == rangetype::pull ){
    AutoSetYRangePull();
  } else if( optype == rangetype::ratio ){
    AutoSetYRangeRatio();
  } else {
  }

  _prevrangetype = optype;
}

/**
 * @brief Y range setting optimized for histogram data presentation.
 *
 * The minimum y axis value is set to be exactly the recorded data minimum
 * (0.3 is TrackY used for all histogram object was set to be automatic).
 * The maximum y axis value is set such that y range of the data only spans
 * the bottom 4/5 of the entire frame, leaving much margin on top for labels
 * and legends. This setting is true for both log and non-log plots, with the
 * expection being that in the case of the log plots, the minimum value is
 * hard set to a minimum of 1e-6.
 */
void
Pad1D::AutoSetYRangeHist()
{
  if( !GetLogy() ){// None log scale
    const double diff = _datamax - _datamin;
    SetYaxisMax( _datamin + diff * 1.25 );
    SetYaxisMin( _datamin );
  } else {
    const double opmin = _datamin < 0 ? 1e-6 : _datamin;
    const double diff  = std::log10( _datamax / opmin );
    SetYaxisMax( opmin * std::pow( 10, diff * 1.25 ) );
    SetYaxisMin( opmin );
  }
}

/**
 * @brief Y axis range setting optimized for scatter data presentation.
 *
 * The data range recorded by the Pad would take up 10/13 of the total
 * y range, 2/13 of the frame would be left for the top margin, and 1/13 would
 * be left for the bottom margin. In case a log plot is requested the "data
 * minimum" would be hard clipped at 1e-6.
 */
void
Pad1D::AutoSetYRangeGraph()
{
  if( !GetLogy() ){// None log scale
    const double diff = _datamax - _datamin;
    SetYaxisMax( _datamin + diff * 1.2 );
    SetYaxisMin( _datamax - diff * 1.1 );
  } else {
    const double opmin = _datamin < 0 ? 1e-6 : _datamin;
    const double diff  = std::log10( _datamax / opmin );
    SetYaxisMax( opmin * std::pow( 10, diff * 1.2 ) );
    SetYaxisMin( _datamax / std::pow( 10, diff*1.1 ) );
    // Minium value is a little smaller
  }
}

/**
 * @brief Y axis range setting optimized for data ratio comparison.
 *
 * The y axis range would be centered around 1. With the range extended to
 * contain the data min/max values, rounded to the closest 0.1 to avoid
 * funny tick values. The maximum range would be 0.4--1.6 to avoid details
 * being lost around data values close to 1. (If there is still nothing in your
 * plot, then there might be something wrong with your data...)
 * An additional 5% margin is left on top and bottom for aesthetic reasons.
 */
void
Pad1D::AutoSetYRangeRatio()
{
  const double diff = std::min( std::max( _datamax - 1, 1-_datamin ), 0.6 );
  // rounding to the closest .1 to avoid weird axis tick labels
  const double opdiff = std::ceil( diff * 10 ) /10.;
  SetYaxisMax( 1 + opdiff *1.05 );
  SetYaxisMin( 1 - opdiff *1.05 );
}

/**
 * @brief Y axis range setting optimized for pull results comparison.
 *
 * The y axis range would be centered around 1. with the range extended to
 * contain the data min/max values, round to the closest 0.1 to avoid funny
 * tick values. The maximum range is hard set to be -5 -- 5 to avoid details
 * around 0 being obscured by extreme values. (If there is nothing in your plot
 * with this range, there might be something wrong with your data...)
 * An additional 10% margin is left on top and bottom for aesthetic reasons.
 */
void
Pad1D::AutoSetYRangePull()
{
  const double diff = std::min( std::max( _datamax, -_datamin ), 5.0 );
  // rounding to the closest .1 to avoid weird axis tick labels
  const double opdiff = std::ceil( diff * 10 ) /10.;
  SetYaxisMax( 1 + opdiff *1.1 );
  SetYaxisMin( 1 - opdiff *1.1 );
}

/**
 * @brief setting axis title according to given unit.
 *
 * Standard format of naming an axis with a certain unit. This fixes the format
 * to "Title [unit]"" (if the unit is not a null string). For the units, both
 * the small and square braces are accepted by the CMS publication standard,
 * the choice of using the square braces is arbitrary.
 */
void
Pad1D::SetAxisTitle(
  TAxis&             axis,
  const std::string& title,
  const std::string& unit
  )
{
  const std::string optitle =
    unit == "" ? title : title + " [" + boost::trim_copy( unit ) +  "]";

  axis.SetTitle( optitle.c_str() );
}

/**
 * @brief Naming x-y axis according to x axis unit and bin-width
 *
 * Method for naming both the X axis and Y axis according to the standard,
 * format. The x axis would be the simples case of "title [unit]", but the y
 * axis would be in the format of "Events/ [bin-width unit]" if the histogram
 * has a uniformed bin width. In the case that the bin-width is not uniformed,
 * then the y axis title should simply be "Events". The final flag allows for
 * forcing the bin-width to be a certain value, and not be calculated
 * automatically.
 *
 * Special not is that the bin-width would at most be set to two digits after
 * the decimal point. If you need more precision, consider scaling your x axis
 * as it makes reading the diagram difficult.
 */
void
Pad1D::SetHistAxisTitles(
  const std::string& title,
  const std::string& unit,
  const double       forcebinwidth )
{
  if( !GetAxisObject() ){
    std::cout << "Error! Axis object doesn't exist yet!" << std::endl;
    return;
  }
  // X axis is straight forwards
  SetAxisTitle( Xaxis(), title, unit );

  // y axis part is more complicated.
  const std::string varwidthfmt = "Events";

  const std::string binwidthstr
    = forcebinwidth != autobinwidth ? forcebinwidth :
      usr::fmt::base::decimal( Xaxis().GetBinWidth( 0 ), -2 );

  // At most 2 digits after decimal point, stripping trailing zeros
  const std::string den
    = binwidthstr == "1" && unit == "" ? "" :
      usr::fstr( "[ %s %s ]", binwidthstr, unit );

  const std::string ytitle
    = Xaxis().IsVariableBinSize()       ? varwidthfmt :
      forcebinwidth == forcevarbinwidth ? varwidthfmt :
      den == ""                         ? varwidthfmt :
      usr::fstr( "Events/ %s", den );

  SetAxisTitle( Yaxis(), ytitle );
}

/**
 * @brief static flag for bin-width determination.
 */
const double Pad1D::autobinwidth = -1;

/**
 * @brief static flag for forcing the variable bin-width scheme.
 */
const double Pad1D::forcevarbinwidth = 0;


}/* plt  */

}/* usr  */
