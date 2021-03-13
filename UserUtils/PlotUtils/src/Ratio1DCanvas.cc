/**
 * @file    Ratio1DCanvas.cc
 * @brief   Implementation of the Ratio1DCanvas functions.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Ratio1DCanvas.hpp"
#else
#include "UserUtils/PlotUtils/Ratio1DCanvas.hpp"
#endif

namespace usr  {

namespace plt {

/**
 * @class PadRatio
 * @details. The pad ratio contains two variable, the ratio between the height
 * of the top and bottom pads; and the gap between the two axis frames relative
 * to the entire canvas.
 */

/**
 * @brief default ratio canvas would have smallist bottom pad, so it is
 * informative but not intrusive.
 */
float PadRatio::default_ratio = 4;

/**
 * @brief default gap would be small such that the two pad having equal x
 * axis scale is immediately obvious, but still look distinct.
 */
float PadRatio::default_gap = 0.005;

/**
 * @{
 * @brief default dimesion for the ratio plots would be a half-page width plot (
 * 0.45 times the full textwidth to be exact), with the top pad remaining
 * having a square frame. And a font size of 8.
 */
length_t Ratio1DCanvas::default_width = 0.45*len::a4textwidth_default();
length_t Ratio1DCanvas::default_height
  = ( 5./4. )*0.45*len::a4textwidth_default();
FontSet Ratio1DCanvas::default_font = FontSet( 8 );
/** @} */

/**
 * @brief Construct a Ratio1DCanvas without specific @ROOT{RooRealVar} to
 * define x axis of both pads.
 */
Ratio1DCanvas::Ratio1DCanvas(
  const length_t  width,
  const length_t  height,
  const PadRatio& pad,
  const FontSet&  set
  ) : Canvas( width, height, set )
{
  Add<Top1DPad>(    PadSize( 0, _splitNDC( pad.ratio ), 1, 1 ) );
  Add<Bottom1DPad>( PadSize( 0, 0, 1, _splitNDC( pad.ratio ) ) );
  _init_margin( pad.gap );
}

/**
 * @brief Construct a Ratio1DCanvas with specific @ROOT{RooRealVar} to define x
 * axis of both pads.
 */
Ratio1DCanvas::Ratio1DCanvas(
  const RangeByVar& range,
  const length_t    width,
  const length_t    height,
  const PadRatio&   pad,
  const FontSet&    set
  ) : Canvas( width, height, set )
{
  Add<Top1DPad>(    PadSize( 0, _splitNDC( pad.ratio ), 1, 1 ), range );
  Add<Bottom1DPad>( PadSize( 0, 0, 1, _splitNDC( pad.ratio ) ), range );
  _init_margin( pad.gap );
  TopPad().SetAxisFont();
  BottomPad().SetAxisFont();
}

Ratio1DCanvas::~Ratio1DCanvas(){}

/**
 * @brief helper function for calculating how to define the PadSize constructor
 * parameter for top and bottom pads.
 */
float
Ratio1DCanvas::_splitNDC( const float x ) const
{
  return ( 1-3.5*Font().lineheight()/Height() )/( x+1. )
         + 3.*Font().lineheight()/Height();
}

/**
 * @brief Initialization the margins of the two pads.
 *
 * @param gap
 */
void
Ratio1DCanvas::_init_margin( const float gap )
{
  SetTopMargin( 1.5*Font().lineheight()/Height() );
  SetLeftMargin( 3.5*Font().lineheight()/Width() );
  SetBottomMargin( 2.0*Font().lineheight()/Height() );
  SetRightMargin(
    std::max(
      ( 3.5*Font().lineheight()/Height() - ( 3.5*Font().lineheight()/Width() ) )
            ,
      ( 0.6*Font().lineheight()/Width() ) )
    );

  // Additional margin setting for gap
  TopPad().SetBottomMargin( 0.5 * gap * Height() / TopPad().AbsHeight() );
  BottomPad().SetTopMargin( 0.5 * gap * Height() / BottomPad().AbsHeight() );
}

/**
 * @{
 * @brief  Margin setting overload to ensure canvas and two pads have common
 * margins for left/right/top/bottom.
 */
void
Ratio1DCanvas::SetTopMargin( const float x )
{
  TopPad().SetTopMargin( x * Height() / TopPad().AbsHeight() );
  TCanvas::SetTopMargin( x );
}

void
Ratio1DCanvas::SetLeftMargin( const float x )
{
  TopPad().SetLeftMargin( x );
  BottomPad().SetLeftMargin( x );
  TCanvas::SetLeftMargin( x );
}

void
Ratio1DCanvas::SetRightMargin( const float x )
{
  TopPad().SetRightMargin( x );
  BottomPad().SetRightMargin( x );
  TCanvas::SetRightMargin( x );
}

void
Ratio1DCanvas::SetBottomMargin( const float x )
{
  BottomPad().SetBottomMargin( x * Height() / BottomPad().AbsHeight() );
  TCanvas::SetBottomMargin( x );
}
/** @} */

/**
 * @brief New font settings for the bottom pad.
 * @details includes hiding the text of the x axis and enlarging the y axis
 * text offset (to account for the smaller default image-size and text size).
 */
void
Top1DPad::SetAxisFont()
{
  Pad1D::SetAxisFont();
  Xaxis().SetTitleSize( 0 );
  Xaxis().SetLabelSize( 0 );
  Yaxis().SetTitleOffset( 2.0 );
}

/**
 * @brief New font settings for the bottom pad.
 *
 * Include large x/y axis text offsets (smaller image size and text size), and
 * changing the Y axis tick marker density for the typical use cases of
 * pull/ratio plots.
 */
void
Bottom1DPad::SetAxisFont()
{
  Pad1D::SetAxisFont();
  Xaxis().SetTitleOffset( 3.8 );
  Yaxis().SetNdivisions( 503 );// Better ticks for axis box
  Yaxis().SetTitleOffset( 2.0 );
}

/**
 * @brief Making the bottom axis object immediately after AxisObject in the top
 * Canvas is found. This is to ensure that the top and bottom bads behave in
 * the same way.
 */
void
Ratio1DCanvas::MakeBottomAxis()
{
  if( BottomPad().GetAxisObject() || !TopPad().GetAxisObject() ){ return; }

  if( TopPad().GetAxisObject()->InheritsFrom( TH1F::Class() ) ){
    TH1F* hist = dynamic_cast<TH1F*>( TopPad().GetAxisObject()->Clone() );
    hist->Reset();
    BottomPad().PlotObj( hist, "AXIS" );
    BottomPad().FrameObj().addObject( hist );
  } else if( TopPad().GetAxisObject()->InheritsFrom( TH1D::Class() ) ){
    TH1F* hist = dynamic_cast<TH1F*>( TopPad().GetAxisObject()->Clone() );
    hist->Reset();
    BottomPad().PlotObj( hist, "AXIS" );
    BottomPad().FrameObj().addObject( hist );
  }

  BottomPad().SetAxisFont();
}

/**
 * @details
 * Given a numerator and denominator histogram, this function generates
 * the scale division results to the two histograms (see static method for
 * details); and claims ownership of the newly generated histogram. The
 * new histogram is then plotted on the bottom pad, with additional plotting
 * method available.
 */
TH1D&
Ratio1DCanvas::PlotScale(
  const TH1D&                   num,
  const TH1D&                   den,
  const std::vector<RooCmdArg>& arglist )
{
  TH1D* ans = ScaleDivide( &num, &den );
  BottomPad().FrameObj().addObject( ans );

  // Parsing arguments.
  const RooArgContainer args( arglist, {
        PlotType( plt::hist ),
        TrackY( TrackY::both )
      } );

  BottomPad().RangeType() = Pad1D::rangetype::ratio;
  BottomPad().PlotHist( ans, args );

  return *ans;
}

/**
 * @details
 * Given a numerator histgram and denominator graph, this function generates the
 * scale division results to the histogram/graph (see static method for details);
 * and claims ownership of the newly generated histogram. The new histogram is
 * then plotted on the bottom pad, with additional plotting method available.
 */
TH1D&
Ratio1DCanvas::PlotScale(
  const TH1D&                   num,
  const TGraph&                 den,
  const std::vector<RooCmdArg>& arglist
  )
{
  // Parsing arguments.
  const RooArgContainer args( arglist, {
        PlotType( plt::hist ),
        TrackY( TrackY::both ),
        ExtrapolateInRatio( kFALSE )
      } );

  // Generating the scaled histogram object
  TH1D* ans = ScaleDivide( &num, &den, 1.0, args.Get<ExtrapolateInRatio>() );

  // Plotting.
  BottomPad().FrameObj().addObject( ans );
  BottomPad().RangeType() = Pad1D::rangetype::ratio;
  BottomPad().PlotHist( ans, args );

  return *ans;
}

/**
 * @details
 * Given a numerator and denominator graph, this function generates
 * the scale division results to the two graphs (see static method for
 * details); and claims ownership of the newly generated graph. The
 * new graph is then plotted on the bottom pad, with additional plotting
 * method available.
 */
TGraphAsymmErrors&
Ratio1DCanvas::PlotScale(
  const TGraph&                 num,
  const TGraph&                 den,
  const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args( arglist, {
        PlotType( plt::simplefunc ),
        TrackY( TrackY::both ),
        ExtrapolateInRatio( kFALSE )
      } );

  TGraphAsymmErrors* ans
    = ScaleDivide( &num, &den, 1.0, args.Get<ExtrapolateInRatio>() );

  BottomPad().FrameObj().addObject( ans );
  BottomPad().RangeType() = Pad1D::rangetype::ratio;
  BottomPad().PlotGraph( ans, args );

  return *ans;
}

/**
 * @details
 * Given a numerator and denominator graph, this function generates the pull
 * results to the two graphs (see static method for details); and claims
 * ownership of the newly generated graph. The new graph is then plotted on the
 * bottom pad, with additional plotting method available.
 */
TGraphAsymmErrors&
Ratio1DCanvas::PlotPull(
  const TGraph&                 num,
  const TGraph&                 den,
  const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args( arglist, {
        PlotType( scatter ),
        TrackY( TrackY::both )
      } );

  TGraphAsymmErrors* ans = PullDivide( &num, &den );

  BottomPad().FrameObj().addObject( ans );
  BottomPad().RangeType() = Pad1D::rangetype::pull;
  BottomPad().PlotGraph( ans, args );

  return *ans;
}

/**
 * @brief dividing two histogram a/b by scaling the numerator by the denominator
 *
 * This function assumes that the numerator and denominator histogram objects
 * have identical binning. The bin content and the bin error of the numerator
 * histogram is then scale to the bin content of the denominator histogram.
 * Bin error of the numerator is ignored.
 *
 * If the bin content of the numerator or the denominator is 0, then the bin
 * content of the resulting histogram would be set to the value cen (1) by
 * default for aesthetic reasons.
 *
 * The output histogram would copy it's styling from the numerator histogram,
 * and the user is responsible for handling the pointer ownership of the
 * generated histogram.
 */
TH1D*
Ratio1DCanvas::ScaleDivide(
  const TH1D*  num,
  const TH1D*  den,
  const double cen
  )
{
  TH1D* ans = new TH1D( *num );

  for( int i = 0; i < num->GetNcells(); ++i ){
    const double n = num->GetBinContent( i );
    const double e = num->GetBinError( i );
    const double d = den->GetBinContent( i );
    if( n == 0 || d == 0 ){
      ans->SetBinContent( i, cen );
      ans->SetBinError( i, 0 );
    } else {
      ans->SetBinContent( i, n/d );
      ans->SetBinError( i, e/d );
    }
  }

  return ans;
}

/**
 * @brief dividing a histogram by a graph. Calculated by scaling the numerator by
 * the denominator
 *
 * This function determines the the denominator by getting the value of the graph
 * at the bin center of the histogram. Uncertainties are also scaled by the same
 * denominator value.
 *
 * If the bin content of the numerator or the denominator is 0, then the bin
 * content of the resulting histogram would be set to the value cen (1) by
 * default for aesthetic reasons.
 *
 * The output histogram would copy it's styling from the numerator histogram, and
 * the user is responsible for handling the pointer ownership of the generated
 * histogram.
 */

TH1D*
Ratio1DCanvas::ScaleDivide(
  const TH1D*   num,
  const TGraph* den,
  const double  cen,
  const bool    extrapolate
  )
{
  TH1D* ans = new TH1D( *num );

  const double xmin = GetXmin( den );
  const double xmax = GetXmax( den );

  for( int i = 0; i < num->GetNcells(); ++i ){
    const double x = num->GetBinCenter( i );
    const double n = num->GetBinContent( i );
    const double e = num->GetBinError( i );
    const double d = den->Eval( x );
    if( n == 0 || d == 0 ){
      ans->SetBinContent( i, cen );
      ans->SetBinError( i, 0 );
    } else if( !extrapolate && ( x < xmin || x > xmax ) ){
      ans->SetBinContent( i, cen );
      ans->SetBinError( i, 0 );
    } else {
      ans->SetBinContent( i, n/d );
      ans->SetBinError( i, e/d );
    }
  }

  return ans;
}

/**
 * @brief Dividing a graph by another graph by scaling the numerator by the
 *        denominator
 *
 * The resulting graph would have the same number of data points as the
 * numerator graph with each of the points y-value scaled by the denominator
 * graph's y-values (linearly interpolated if the two graph's x values doesn't
 * exactly match). The y-errors of the numerator graph will also be scaled but
 * the x-errors and x-value would be left unchanged. The error bars of the
 * denominator graph would be ignored.
 *
 * The output graph would copy it's styling from the numerator graph, and the
 * user is responsible for handling the pointer ownership of the generated
 * graph.
 */
TGraphAsymmErrors*
Ratio1DCanvas::ScaleDivide(
  const TGraph* num,
  const TGraph* den,
  const double  cen,
  const bool    extrapolate  )
{
  const double xmin = GetXmin( den );
  const double xmax = GetXmax( den );
  std::vector<double> x_list;
  std::vector<double> y_list;
  std::vector<double> xerrlo_list;
  std::vector<double> xerrhi_list;
  std::vector<double> yerrlo_list;
  std::vector<double> yerrhi_list;

  for( int i = 0; i < num->GetN(); ++i ){
    const double origx  = num->GetX()[i];
    const double origy  = num->GetY()[i];
    const double xerrlo = num->GetErrorXlow( i );
    const double xerrhi = num->GetErrorXhigh( i );
    const double yerrlo = num->GetErrorYlow( i );
    const double yerrhi = num->GetErrorYhigh( i );

    const double deny = den->Eval( origx );

    if( !extrapolate && ( origx < xmin || origx > xmax ) ){
      // Do nothing
    } else if( deny == 0 || origy == 0 ){
      x_list.push_back( origx );
      y_list.push_back( cen );
      xerrlo_list.push_back( 0 );
      xerrhi_list.push_back( 0 );
      yerrlo_list.push_back( 0 );
      yerrhi_list.push_back( 0 );
    } else {
      x_list.push_back( origx );
      y_list.push_back( origy/deny );
      xerrlo_list.push_back( xerrlo );
      xerrhi_list.push_back( xerrhi );
      yerrlo_list.push_back( yerrlo/deny );
      yerrhi_list.push_back( yerrhi/deny );
    }

    if( y_list.back() < 0.8 ){
      std::cout << i << " "
                << origx << " "
                << origy << std::endl;
    }
  }

  assert( x_list.size() == y_list.size() );
  assert( x_list.size() == xerrlo_list.size() );
  assert( x_list.size() == xerrhi_list.size() );
  assert( x_list.size() == yerrlo_list.size() );
  assert( x_list.size() == yerrhi_list.size() );

  TGraphAsymmErrors* ans
    = new TGraphAsymmErrors( x_list.size()
                           , x_list.data(), y_list.data()
                           , xerrlo_list.data(), xerrhi_list.data()
                           , yerrlo_list.data(), yerrhi_list.data() );

  ans->SetTitle( "" );
  ans->GetXaxis()->SetLimits(
    num->GetXaxis()->GetXmin(),
    num->GetXaxis()->GetXmax() );
  // duplicating style
  ans->SetLineColor( num->GetLineColor() );
  ans->SetLineStyle( num->GetLineStyle() );
  ans->SetLineWidth( num->GetLineWidth() );
  ans->SetFillColor( num->GetFillColor() );
  ans->SetFillStyle( num->GetFillStyle() );
  ans->SetMarkerColor( num->GetMarkerColor() );
  ans->SetMarkerStyle( num->GetMarkerStyle() );
  ans->SetMarkerSize( num->GetMarkerSize() );

  return ans;
}

/**
 * @brief Generating the pull diagram when comparing two graphs
 *
 * The resulting graph would have the following value:
 * - gen.Y(x) = num.Y(x) - den.Y(x)
 * - gen.ErrY(x) = num.ErrY(x) / den.ErrY(x)
 * - gen.ErrX(x) = 0
 *
 * Interpolation is used on the denominator to obtain the error and value.
 *
 * In the special case that den.ErrY(x) is zero, gen.Y(x) is set to the cen value
 * (defaults to 0 for aesthetic reasons).
 */
TGraphAsymmErrors*
Ratio1DCanvas::PullDivide(
  const TGraph* num,
  const TGraph* den,
  const double  cen )
{
  TGraphAsymmErrors* ans = new TGraphAsymmErrors( num->GetN() );

  // duplicating style
  ans->SetLineColor( num->GetLineColor() );
  ans->SetLineStyle( num->GetLineStyle() );
  ans->SetLineWidth( num->GetLineWidth() );
  ans->SetFillColor( num->GetFillColor() );
  ans->SetFillStyle( num->GetFillStyle() );
  ans->SetMarkerColor( num->GetMarkerColor() );
  ans->SetMarkerStyle( num->GetMarkerStyle() );
  ans->SetMarkerSize( num->GetMarkerSize() );

  // Making error graphs for interpolation
  TGraph denerrup   = TGraph( den->GetN() );
  TGraph denerrdown = TGraph( den->GetN() );

  for( int i = 0; i < den->GetN(); ++i ){
    denerrup.SetPoint( i, den->GetX()[i], den->GetErrorYhigh( i ) );
    denerrdown.SetPoint( i, den->GetX()[i], den->GetErrorYlow( i ) );
  }

  for( int i = 0; i < num->GetN(); ++i ){
    const double x    = num->GetX()[i];
    const double diff = num->GetY()[i] - den->Eval( x );
    const double err  = diff > 0 ? denerrup.Eval( x ) :
                        denerrdown.Eval( x );

    ans->SetPointEXhigh( i, 0 );
    ans->SetPointEXlow( i, 0 );
    if( err > 0  ){
      ans->SetPoint( i, x, diff / err );
      ans->SetPointEYhigh( i, num->GetErrorYhigh( i ) / err );
      ans->SetPointEYlow( i, num->GetErrorYlow( i ) / err );
    } else {
      ans->SetPoint( i, x, cen );
      ans->SetPointEYhigh( i, 0 );
      ans->SetPointEYlow( i, 0 );
    }
  }


  return ans;
}

}/* plt */

}/* usr  */
