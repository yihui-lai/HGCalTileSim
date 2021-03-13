/**
 * @file    Simple1DCanvas.cc
 * @brief   Implementation of the Simple1DCanvas class.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Simple1DCanvas.hpp"
#else
#include "UserUtils/PlotUtils/Simple1DCanvas.hpp"
#endif

namespace usr  {

namespace plt  {

/**
 * @{
 * @brief default construction dimensions
 *
 * A square of 0.8 times the text box width in standard CMS documents. and the
 * font size of 11pt (the same as the main text size used in CMS documents).
 */
length_t Simple1DCanvas::default_width  = 0.8 * len::a4textwidth_default();
length_t Simple1DCanvas::default_height = 0.8 * len::a4textwidth_default();
FontSet  Simple1DCanvas::default_font = FontSet();
/** @} */

/**
 * @brief Construct a Simple1DCanvas without a specific @ROOT{RooRealVar} to
 * define the x axis.
 */
Simple1DCanvas::Simple1DCanvas(
  const length_t width,
  const length_t height,
  const FontSet& fontset
  ) : Canvas( width, height, fontset )
{
  Add<Pad1D>( PadSize( 0, 0, 1, 1 ) );
  _init_margin();
}

/**
 * @brief Constructing a Simple1DCanvas with a specific @ROOT{RooRealVar} to
 * define the x axis.
 */
Simple1DCanvas::Simple1DCanvas(
  const RangeByVar& range,
  const length_t    width,
  const length_t    height,
  const FontSet&    fontset
  ) : Canvas( width, height, fontset )
{
  Add<Pad1D>( PadSize( 0, 0, 1, 1 ), range );
  _init_margin();
}

Simple1DCanvas::~Simple1DCanvas(){}

/**
 * First defining the default margins on the object:
 * - 3.5 line-heights for the left margin
 * - 2.5 line-heights for the bottom margin.
 * - 1.5 lineheight for the top margin.
 * - The right margin is automatically assigned such that the framebox would be
 *   the same ratio as the canvas. If such a margin value would be greater than
 *   0.3 lineheights, the margin would be set to 0.3 line-heights to avoid
 *   awkwardly small margins.
 */
void
Simple1DCanvas::_init_margin()
{
  Pad().SetTopMargin(    1.5*Font().lineheight()/Height() );
  Pad().SetLeftMargin(   3.5*Font().lineheight()/Width() );
  Pad().SetBottomMargin( 2.5*Font().lineheight()/Height() );
  Pad().SetRightMargin( std::max(
    ( 4.0*Font().lineheight()/Height() - ( 3.5*Font().lineheight()/Width() ) ),
    ( 0.3*Font().lineheight()/Width() ) )
    );
}

}/* plt  */

}/* usr  */
