/**
 * @file    Pad1D_Draw.cc
 * @brief   Implementing the drawing functions that spawn new objects.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Pad1D.hpp"
#else
#include "UserUtils/PlotUtils/Pad1D.hpp"
#endif

#include "CmdSetAttr.hpp"
#include <boost/format.hpp>


namespace usr  {

namespace plt  {

/**
 * @{
 * @brief Function for getting a nice text placement in the corner of a Pad.
 *
 * An extra margin based on the current axis tick-length would be added.
 */
float
Pad1D::InnerTextLeft() const
{
  return GetLeftMargin() + 1.5 * Yaxis().GetTickLength();
}

float
Pad1D::InnerTextRight() const
{
  return 1. - GetRightMargin() - 1.5 * Yaxis().GetTickLength();
}

float
Pad1D::InnerTextHCenter() const
{
  return ( InnerTextLeft() + InnerTextRight() ) /2;
}

float
Pad1D::InnerTextTop() const
{
  return 1 - GetTopMargin() - 1.5 * Xaxis().GetTickLength();
}

float
Pad1D::InnerTextBottom() const
{
  return GetBottomMargin() + 1.5 * Xaxis().GetTickLength();
}

float
Pad1D::InnerTextVCenter() const
{
  return ( InnerTextTop() + InnerTextBottom() ) /2;
}
/** @} */

/**
 * @brief Drawing horizontal line
 *
 * Drawing a horizontal line across the pad at a fixed y value. Accepts the
 * LineColor, LineStyle, LineWidth options for change the styling of the line.
 * Also returns a reference to the line generated.
 */
TLine&
Pad1D::DrawHLine(
  const double&                 y,
  const std::vector<RooCmdArg>& arglist  )
{
  RooArgContainer args( arglist );
  TLine& line = MakeObj<TLine>( GetXaxisMin(), y, GetXaxisMax(), y );
  SetLineAttr( line, args );
  PadBase::PlotObj( line );
  return line;
}

/**
 * @brief Drawing vertical line
 *
 * Drawing a vertical line across the pad at a fixed x value. Accepts the
 * LineColor, LineStyle, LineWidth options for change the styling of the line.
 * Also return a reference to the line generated.
 */
TLine&
Pad1D::DrawVLine(
  const double&                 x,
  const std::vector<RooCmdArg>& arglist )
{
  RooArgContainer args( arglist );
  TLine& line = MakeObj<TLine>(
    x, Yaxis().GetXmin(),
    x, Yaxis().GetXmax() );
  SetLineAttr( line, args );
  PadBase::PlotObj( line );
  _linelist.push_back( &line );
  return line;
}

/**
 * @brief Drawing CMS label
 *
 * Drawing the CMS label inside the frame at the top left corner.
 * Additional styling will include:
 * - The font size of the label will be a "large" variant to the Canvas settings
 * - The Top left corner of the text would be set to 1.5 times the tick length.
 * - The man label (CMS) would be set to bold
 * - The second label (if defined, such as "Preliminary") would be set to
 *   italic.
 * @param tag  The secondary tag to be displayed below the CMS label font.
 * @param main The main tag to be plotted (default to CMS but could be changed
 *             if you are working for another experiment)
 */
void
Pad1D::DrawCMSLabel( const std::string& tag, const std::string& main )
{
  SetTextCursor( InnerTextLeft(), InnerTextTop(), font::top_left );
  WriteLine( "#bf{"+main+"}", TextSize( Font().large() ) );
  if( tag != "" ){
    WriteLine( "#it{" + tag + "}", TextSize( Font().large() ) );
  }
}

/**
 * @brief Drawing luminosity label.
 *
 * Writing the Luminosity label outside of the frame on the top-right corner.
 * This string in directly pass to the TLatex object, so all styling here would
 * be user defined. The font setting would be identical to that of the canvas
 * setting.
 * @param x the string to use as the luminosity label
 */
void
Pad1D::DrawLuminosity( const std::string& x )
{
  const float ymargin = 0.3*Xaxis().GetTickLength();
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
 * @param lumi The luminosity in units of fb-1
 */
void
Pad1D::DrawLuminosity( const double lumi )
{
  DrawLuminosity( usr::fstr( "%.1f fb^{-1} (13 TeV)", lumi ) );
}

/**
 * @brief Resetting the vertical lines according to the data plotted on the pad
 *
 */
void
Pad1D::FixVLines()
{
  for( auto line : _linelist ){
    line->SetY1( GetYaxisMin() );
    line->SetY2( GetYaxisMax() );
  }
}

}/* plt  */

}/* usr  */
