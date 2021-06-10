#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Flat2DCanvas.hpp"
#else
#include "UserUtils/PlotUtils/Flat2DCanvas.hpp"
#endif

namespace usr {

namespace plt {

/**
 * @{
 * @brief default construction dimensions
 *
 * The height of the canvas is 0.5 times the width  of a standard CMS, the width
 * is 1.6 times the height. font size of 11pt (the same as the main text size
 * used in CMS documents). The margin will be set up so that the box of the plot
 * is square.
 **/
length_t Flat2DCanvas::default_width  = 1.0  * len::a4textwidth_default();
length_t Flat2DCanvas::default_height = 0.65 * len::a4textwidth_default();
FontSet Flat2DCanvas::default_font    = FontSet();
/** @} */


/**
 * @brief Construct a Flat2DCanvas.
 */
Flat2DCanvas::Flat2DCanvas(
  const length_t width,
  const length_t height,
  const FontSet& fontset ) :
  Canvas( width, height, fontset )
{
  Add<Pad2DFlat>( PadSize( 0, 0, 1, 1 ) );
  init_margin();
}

Flat2DCanvas::~Flat2DCanvas(){}

/**
 * First defining the default margins on the object:
 * - 2.5 line-heights for the top, bottom, and right margin.
 * - Left margin is set so the plot box remains a square. (Pad2DFlat is designed
 *   so that legend will be drawn on the left)
 */
void
Flat2DCanvas::init_margin()
{
  Pad().SetTopMargin(    1.5*Font().lineheight()/Height() );
  Pad().SetRightMargin(  5.5*Font().lineheight()/Width() );
  Pad().SetBottomMargin( 2.5*Font().lineheight()/Height() );
  Pad().SetLeftMargin(
    1 - double(Height())/Width() - 1.5*Font().lineheight()/Width()
    );
}

}

}
