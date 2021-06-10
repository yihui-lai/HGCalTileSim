#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Pad2DFlat.hpp"
#else
#include "UserUtils/PlotUtils/Pad2DFlat.hpp"
#endif


#include <boost/format.hpp>

namespace usr {

namespace plt {

/**
 * @brief Drawing CMS label
 *
 * Drawing the CMS label inside the frame outside the top left corner.
 * Additional styling will include:
 * - The font size of the label will be a "large" variant to the Canvas settings
 * - The bottom left corner would be aligned to the top left corner of the frame
 *   corner, with a vertical offset of 0.25 tick-lengths
 * - The main label (CMS by default) would be set to bold
 * - The second label (if defined, such as "Preliminary") would be set to
 *   italic.
 * @param tag  The secondary tag to be displayed beside the CMS label font.
 * @param main The main tag to be plotted (default to CMS but could be changed
 *             if you are working for another experiment)
 */
void
Pad2DFlat::DrawCMSLabel( const std::string& tag, const std::string& main )
{
  const float tmp     = FontSize();
  const float ymargin = 0.25 * Xaxis().GetTickLength();
  const float cursorx = GetLeftMargin();
  const float cursory = 1 - GetTopMargin() + ymargin;
  _latex.SetTextSize( Font().large() );

  SetTextCursor( cursorx, cursory, font::bottom_left );
  WriteLine( usr::fstr("#bf{%s} #it{%s}", main, tag ) );
  _latex.SetTextSize( tmp );
}


/**
 * @brief Drawing luminosity label.
 *
 * Writing the Luminosity label outside of the frame on the top-right corner.
 * This string in directly pass to the TLatex object, so all styling here would
 * be user defined. The font setting would be identical to that of the canvas
 * setting.
 */
void
Pad2DFlat::DrawLuminosity( const std::string& x )
{
  const float ymargin = 0.25 * Xaxis().GetTickLength();
  const float cursorx = 1 - GetRightMargin();
  const float cursory = 1 - GetTopMargin() + ymargin;

  SetTextCursor( cursorx, cursory, font::bottom_right );
  WriteLine( x );
}

/**
 * @brief Drawing luminosity label.
 *
 * Writing the luminosity label in the standard format of
 * "XX.X fb^{-1} (13TeV)"
 */
void
Pad2DFlat::DrawLuminosity( const double lumi )
{
  DrawLuminosity( usr::fstr( "%.1f fb^{-1} (13 TeV)", lumi ) );
}

}

}
