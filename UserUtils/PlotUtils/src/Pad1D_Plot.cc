/**
 * @file    Pad1D_Plot.cc
 * @brief   Implementing the data plotting functions.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/Maths.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/interface/Miscellaneous.hpp"
#include "UserUtils/PlotUtils/interface/Pad1D.hpp"
#else
#include "UserUtils/Common/Maths.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/MathUtils/Miscellaneous.hpp"
#include "UserUtils/PlotUtils/Pad1D.hpp"
#endif

#include <limits>
#include <random>

#include "CmdSetAttr.hpp"

#include "TDecompChol.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TList.h"

// static variables for new object generation
static const std::string genaxisname = "axishist";

namespace usr  {

namespace plt {

/**
 * Plotting 1D histograms accepts these options:
 * - PlotType: Defining how the binned data should be represented on the pad,
 *   following types are supported
 *   - plottype::hist (default) - standard square curve/block diagram for
 *     displaying a diagram, equivalent to the `"HIST"` options provided in the
 *     TH1 object.
 *   - plottype::scatter - Plotting points with error bars. Complying to the CMS
 *     plotting convention with the horizontal error bar being suppressed for
 *     fixed bin-width data.
 *   - plottype::histstack - Stacking this histogram into the histogram stack on
 *     the Pad. If the stack doesn't already exist, a new stack is created. Note
 *     that all histograms in the stack will be plotted with the hist style and
 *     cannot be changed. Notice that the stack must be plotted consecutively,
 *     meaning that any non histstack options would cause the stack to be
 *     finalized and plotted onto the Pad.
 *   - plottype::histnewstack - Force the creation of the new histogram stack,
 *     used only if you are plotting two histogram stacks in the same Pad.
 *   - plottype::histerr - Drawing histogram uncertainty as a shaded box region.
 *   - Raw string: The user can use the ROOT style string options to define how
 *     the Histogram options should be plotted. The strings `"SAME"` and `"HIST"`
 *     would be handed by the functions and will be stripped.
 *
 * - EntryText: String to add in the legend entry. If this options is not
 *   present, then this object will NOT appear in the generated legend. Notice
 *   that the attributes to display in the legend would be generated from the
 *   PlotType used.
 *
 * - TrackY: Whether or not the y-axis range should be adjusted according to the
 *   newly added histogram. By default, only the max values of the histogram
 *   would be used.
 *
 * - PlotUnder: Specifying the that plot object be plotted underneath a specific
 *   object. Not that if the object specified in PlotUnder is not found in the
 *   pad, this option will have no effect and will raise no excpetions.
 *
 * One side note is that fitted functions will have its DrawOptions cleared from
 * the histogram! The user should be the one explicitly invoking the plotting
 * behavior, given more freedom to what the final plot will look like (See the
 * Pad1D::PlotFunc() method for more details). This design aspect is in stark
 * contrast with the design of ROOT objects. So beware of unwanted behaviour.
 */
TH1D&
Pad1D::PlotHist( TH1D& obj, const std::vector<RooCmdArg>& arglist )
{
  // If TProfile is given, send to PlotProfile function instead
  if( obj.InheritsFrom( TProfile::Class() ) ){
    return PlotProfile( dynamic_cast<TProfile&>( obj ), arglist );
  }

  // Getting the flags
  const RooArgContainer args( arglist, {
        PlotType( hist ),
        TrackY( TrackY::max ),
        FillColor( 0, 0 )
      } );

  if( !GetAxisObject() ){
    // MUST ust a clone, otherwise messes with THStack
    TH1D& axisobj = MakeObj<TH1D>( obj );
    axisobj.Reset();
    axisobj.SetStats( 0 );
    axisobj.SetTitle( "" );
    axisobj.SetName( ( genaxisname+RandomString( 6 ) ).c_str() );
    PlotObj( axisobj, "AXIS" );
    this->SetAxisFont();
  }

  // Forcing no statistics. and wiping title
  obj.SetStats( 0 );
  obj.SetTitle( "" );

  // Removing the poly marker from TSpectrum search function
  TObject* polymarker = obj.GetListOfFunctions()->FindObject( "TPolyMarker" );
  obj.GetListOfFunctions()->RecursiveRemove( polymarker );

  // Forcing fit functions to not be drawn
  for( const auto&& func : *( obj.GetListOfFunctions() ) ){
    func->SetBit( TF1::kNotDraw, true );
  }

  // Running the draw commands.
  const int pt = args.Get<PlotType>();

  // Flushing the _working stack if hist is no longer used
  if( pt != histstack && _workingstack ){
    PlotObj( _workingstack, "HIST SAME NOCLEAR" );
    _workingstack = nullptr;
  }

  // Plot case for different plot types
  switch( pt ){
  case plottype::hist:
    PlotObj( obj, "HIST SAME" );
    break;
  case plottype::scatter:
    obj.GetXaxis()->IsVariableBinSize() ?
    PlotObj( obj,  "P L E SAME" ) :
    PlotObj( obj, "P E X0 SAME" );
    break;
  case plottype::histerr:
    PlotObj( obj,     "E2 SAME" );
    break;
  case plottype::histstack:
    if( _workingstack  == 0 ){
      _workingstack = &MakeObj<THStack>(
        ( "stack" + RandomString( 12 ) ).c_str(), "" );
    }
    _workingstack->Add( &obj, "HIST" );
    break;
  case plottype::histnewstack:
    _workingstack = &MakeObj<THStack>(
      ( "stack" + RandomString( 12 ) ).c_str(), "" );
    _workingstack->Add( &obj, "HIST" );
    break;
  case plottype::plottype_dummy:
    PlotObj( obj, ( args.Get<PlotType>().str()+" SAME" ).c_str() );
    break;
  default:
    std::cerr << "Skipping over invalid value ("<< pt <<")" << std::endl;
  }

  if( _workingstack ){
    TrackObjectY( *_workingstack, args.Get<TrackY>() );
  } else {
    TrackObjectY( obj, args.Get<TrackY>() );
  }

  // Adding legend
  if( args.Has<EntryText>() ){
    AddLegendEntry( obj, args.Get<EntryText>(), args.Get<PlotType>() );
  }

  // Moving to under something
  if( args.Has<PlotUnder>() ){
    if( _workingstack ){
      PadBase::MoveTargetToBefore( *_workingstack, args.Get<PlotUnder>() );
    } else {
      PadBase::MoveTargetToBefore( obj, args.Get<PlotUnder>() );
    }
  }

  SetLineAttr( obj, args );
  SetFillAttr( obj, args );
  SetMarkAttr( obj, args );

  return obj;
}

/**
 * @brief  Plotting a TProfile object
 *
 * Due to inheritance between TProfile and the TH1D, the original PlotHist
 * doesn't work with TProfile. Instead, this functions generates a new TH1D from
 * the TProfile using the ProjectionX function. The return value would then be
 * the newly generated histogram. So beware of unexpected behaviour.
 */
TH1D&
Pad1D::PlotProfile( TProfile& obj, const std::vector<RooCmdArg>& args )
{
  // Generating new profile objects
  TH1D* _newhist = obj.ProjectionX( usr::RandomString( 6 ).c_str(), "E" );
  ClaimObject( _newhist );
  return PlotHist( _newhist, args );
}

/**
 * Plotting of the TGraph object has the following supporting options:
 *
 * - `PlotType`: Defining how the data should be represented on the Pad. The
 *   supported types are:
 *   - `plottype::simplefunc`(default): single polyline joining the data points.
 *     This is used if the graph represents a function sample.
 *   - `plottype::fittedfunc`: a polyline joining the data points with the Y
 *     error represented with a shaded region. This is used if the graph
 *     represents a function sample with additional sampling for fitting
 *     uncertainties
 *   - `plottype::scatter`: plotting the data as data points with error bars.
 *     Unlike the histograms, this would not attempt to adjust the x error bars.
 *
 * - EntryText: String to add in the legend entry. If this options is not
 *   present, then this object will NOT appear in the generated legend. Notice
 *   that the attributes to display in the legend would be generated from the
 *   PlotType used.
 *
 * - TrackY: Whether or not the y-axis range should be adjusted according to the
 *   newly added graph. By default, both the maximum and minimum will be tracked
 *   if this graph is the first thing to be plotted on the pad (excluding the
 *   histogram used for handling the axis), otherwise, nothing will be tracked.
 *
 * - PlotUnder: Specifying the that plot object be plotted underneath a specific
 *   object. Not that if the object specified in PlotUnder is not found in the
 *   pad, this option will have no effect and will raise no excpetions.
 *
 * One side note is that fitted functions will have its DrawOptions cleared from
 * the histogram! The user should be the one explicitly invoking the plotting
 */
TGraph&
Pad1D::PlotGraph( TGraph& obj, const std::vector<RooCmdArg>& arglist )
{
  // Early Exit for Graphs without any data points
  if( obj.GetN() <= 0 ){
    std::cerr << "Cannot plot TGraphs with no data points!" << std::endl;
    return obj;
  }

  // Getting the flags
  const RooArgContainer args( arglist,
      {
        PlotType( simplefunc ),
        !GetAxisObject() ? TrackY( TrackY::both ) : TrackY( TrackY::none ),
        FillColor( 0, 0 )
      } );

  // If no axis are available. Generating a TH1 object for axis:
  if( !GetAxisObject() ){
    auto& axishist = MakeObj<TH1D>(
      ( genaxisname + RandomString( 6 ) ).c_str(),
      "",
      10, GetXmin( obj ), GetXmax( obj ) );
    axishist.SetStats( 0 );
    PadBase::PlotObj( axishist, "AXIS" );
    SetAxisFont();
  }

  // Object fixing
  obj.SetTitle( "" );// Forcing clear title. This should be handled by Canvas.

  // Forcing fit functions to not be drawn
  for( const auto&& func : *( obj.GetListOfFunctions() ) ){
    func->SetBit( TF1::kNotDraw, true );
  }

  const int pt = args.Get<PlotType>();

  switch( pt ){
  case simplefunc:
    PadBase::PlotObj( obj, "LX" );
    break;
  case fittedfunc:
    PadBase::PlotObj( obj,  "3" );// Draw Error with fill region and then
    PadBase::PlotObj( obj, "LX" );// Draw the central line. All errors disabled.
    break;
  case scatter:
    // Points, no error bar end ticks, show error bar for points  outside range.
    PadBase::PlotObj( obj, "PZ0" );
    break;
  case plottype_dummy:
    PlotObj( obj, args.Get<PlotType>().c_str() );
    break;

  default:
    std::cerr << "Skipping over invalid value" << std::endl;
    break;
  }

  TrackObjectY( obj, args.Get<TrackY>() );

  // Adding legend
  if( args.Has<EntryText>() ){
    AddLegendEntry( obj, args.Get<EntryText>(), args.Get<PlotType>() );
  }

  // Moving to under something
  if( args.Has<PlotUnder>() ){
    PadBase::MoveTargetToBefore( obj, args.Get<PlotUnder>() );
  }

  // Setting styling attributes
  SetLineAttr( obj, args );
  SetFillAttr( obj, args );
  SetMarkAttr( obj, args );

  return obj;
}

/**
 * Plotting a TF1 object is done by generating a TGraph with 300 samples points
 * across the x axis, and plotting the TGraph instead. All TGraph plotting
 * options will be available for the TF1 plotting. There is a new plotting
 * options VisualizeError, which generates a TGraphErrors graph by randomly
 * sampling the parameter space according to the correlation matrix given by a
 * @ROOT{TFitResult}.
 */
TGraph&
Pad1D::PlotFunc( TF1& func, const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args(
    arglist,
      {
        PlotType(// Inserting default plotting style
          RooArgContainer::CheckList( arglist, VisualizeError::CmdName ) ?
          fittedfunc : simplefunc ),
        RooFit::Precision( 1e-3 )
      }
    );

  // If no axis are available. Generating a TH1 object for axis:
  if( !GetAxisObject() ){
    auto& axishist = MakeObj<TH1D>(
      ( genaxisname + RandomString( 10 ) ).c_str(),
      "",
      10, func.GetXmin(), func.GetXmax() );
    axishist.SetStats( 0 );
    PadBase::PlotObj( axishist, "AXIS" );
    SetAxisFont();
  }

  return PlotGraph( MakeTF1Graph( func, args ), args );
}

/**
 * @brief Generating a graph representation of the TF1 object.
 */
TGraph&
Pad1D::MakeTF1Graph( TF1& func, const RooArgContainer& args  )
{
  static const unsigned psample = 300;
  const std::string graphname   = func.GetName()
                                  + std::string( "_gengraph" )
                                  + RandomString( 6 );

  const double xmax      = func.GetXmax();
  const double xmin      = func.GetXmin();
  const double xspace    = args.Get( "Precision" ).getDouble( 0 );
  const unsigned xsample = 1 / ( xspace ) + 1;

  std::vector<double> x( xsample );
  std::vector<double> y( xsample );
  std::vector<double> yerrhi( xsample );
  std::vector<double> yerrlo( xsample );
  std::vector<double> zero( xsample );

  // Getting common elements for graph generation
  for( unsigned i = 0; i < xsample; ++i ){
    const double xval = xmin + i * ( xmax-xmin ) * xspace;
    const double yval = func.Eval( xval );

    x[i]      = xval;
    y[i]      = yval;
    yerrlo[i] = 0;
    yerrhi[i] = 0;
    zero[i]   = 0;
  }

  if( !args.Has<VisualizeError>() ){
    TGraph& graph = MakeObj<TGraph>( x.size(), x.data(), y.data() );
    graph.SetName( graphname.c_str() );
    return graph;
  } else {
    const auto cmd    = args.Get<VisualizeError>();
    const auto& fit   = cmd.GetTFitResult();
    const double zval = cmd.getDouble( 0 );

    const std::vector<double> bestfit_param = fit.Parameters();

    // Getting matrix for random parameter generation
    const TMatrixDSym cormatrix = fit.GetCovarianceMatrix();
    const TMatrixD tmatrix      = usr::DecompCorvariance( cormatrix );

    TVectorD vec( tmatrix.GetNrows() );
    std::mt19937 gen;
    std::normal_distribution pdf( 0.0, 1.0 );

    // Random sample for parameter space
    for( unsigned i = 0; i < psample; ++i ){

      // Generating random variation using gaussian
      for( int j = 0; j < vec.GetNrows(); ++j ){
        vec[j] = pdf( gen );
      }

      // Forcing the vector onto unit sphere, then transforming according to
      // the covariance matrix
      vec = ( zval /sqrt( vec.Norm2Sqr() ) ) * vec;
      vec = tmatrix * vec;

      // Shifting to central value of function.
      for( int j = 0; j < vec.GetNrows(); ++j ){
        func.SetParameter( j, vec[j] + bestfit_param[j] );
      }

      // Finding evelope of randomly generated parameter values
      for( unsigned j = 0; j < xsample; ++j ){
        const double xval = x.at( j );
        const double yerr = func.Eval( xval ) - y.at( j );
        yerrhi[j] = std::max( yerr, yerrhi.at( j ) );
        yerrlo[j] = std::max( -yerr, yerrlo.at( j ) );
      }
    }

    // Reseting the function parameters to best value
    func.SetParameters( bestfit_param.data() );

    TGraphAsymmErrors& graph = MakeObj<TGraphAsymmErrors>( x.size(),
      x.data(), y.data(),
      zero.data(), zero.data(),
      yerrlo.data(), yerrhi.data() );

    graph.SetName( graphname.c_str() );
    return graph;
  }
}

/**
 * Function for plotting all of the RooAbsData objects, the supported options
 * allowed include:
 * - *Any* RooCmdArg supported by the RooAbsData::plotOn function. These will
 *   take precedence to custom defined options. With the excpetion of
 *   RooFit::DrawOptions.
 *
 * - PlotType: Specifying the representation for the data object on the Pad. Any
 *   of the types valid for the PlotGraph functions would be valid. By default,
 *   the scatter plot method would be used.
 *
 * - EntryText: String to add in the legend entry. If this options is not
 *   present, then this object will NOT appear in the generated legend. Notice
 *   that the attributes to display in the legend would be generated from the
 *   PlotType used. Notice that this method is overwritten by the
 *   RooFit::Invisible option, which forces the object to be hidden.
 *
 * - TrackY: Whether or not the y-axis range should be adjusted according to the
 *   newly added graph. By default, only the maximum value would be used.
 *
 * Additional automation of options include suppressing the X error bars of the
 * generated graph if an equal binning is used for the data set and the
 * RooFit::XErrorSize is not specified. This is to align the plotting style
 * conventions within the CMS collaboration. Another change to the graphs is that
 * the y uncertainties in the zero bins are suppressed, again, to align with the
 * plotting style conventions.
 *
 * The function returns a reference to the generated TGraph object by the
 * underlying RooPlot object, The object would be owned by the RooPlot and be
 * destroyed when the Pad1D goes out of scope.
 */
TGraphAsymmErrors&
Pad1D::PlotData(
  RooAbsData&                   data,
  const std::vector<RooCmdArg>& arglist
  )
{
  const RooArgContainer args(
    arglist,
      {
        PlotType( scatter ),
        TrackY( TrackY::max )
      }
    );

  TGraphAsymmErrors& ans = MakeDataGraph( data, args );

  if( !args.Has( "Invisible" ) ){
    if( _prevrangetype == rangetype::aut ){
      // Forcing the range type to be histogram like.
      _prevrangetype = rangetype::hist;
    }
    PlotGraph( ans, args );
  }

  return ans;
}

/**
 * @brief Helper method for generating the a TGraph representing a RooAbsData
 * plot.
 *
 * This will detect the binning method specified by the various options (using
 * the present binning of the RooFrame variable otherwise), and in the case of
 * uniform binning, suppress the x axis error bars, as specified by the CMS
 * publication conventions. Additionally, Bins with zero entires would have their
 * y error bars suppressed.
 */
TGraphAsymmErrors&
Pad1D::MakeDataGraph( RooAbsData&            data,
                      const RooArgContainer& args )
{
  // RooCmdArg for suppressing the X errors in the final graph
  static const RooCmdArg suppressxerror = RooFit::XErrorSize( 0 );
  // Generating the requirements for actual RooFit PlotOn calls
  RooLinkedList oplist = args.MakeRooList();

  auto IsUniform
    = [this]( const RooCmdArg& cmd ) -> bool {
        // Double assignment will always be uniform
        if( cmd.getDouble( 0 ) != cmd.getDouble( 1 ) ){ return true;}
        // Getting plot variable bining scheme
        if( cmd.getString( 0 ) ){
          const auto binname = cmd.getString( 0 );
          const auto plotvar = this->_frame.getPlotVar();
          if( plotvar->hasBinning( binname ) ){
            return plotvar->getBinning( binname ).isUniform();
          }
        } else if( cmd.getObject( 0 ) ){
          auto binobj = dynamic_cast<const RooAbsBinning*>( cmd.getObject( 0 ) );
          if( binobj  ){
            return binobj->isUniform();
          }
        }
        return true;// Returning true by default
      };

  // Option for suppressing x error bars
  if( !args.Has( "Binning" ) ){
    if( _frame.getPlotVar()->getBinning().isUniform() ){
      oplist.Add( suppressxerror.Clone() );
    }
  } else if( IsUniform( args.Get( "Binning" ) ) ){
    oplist.Add( suppressxerror.Clone() );
  }

  // Generating the object via RooFit
  TGraphAsymmErrors& ans = GenGraph( data, oplist );

  // Additional Fixes for RooFit generated objects
  for( int i = 0; i < ans.GetN(); ++i ){
    if( ans.GetY()[i] == 0 ){// Suppressing errors in zero bins
      ans.SetPointError( i, 0, 0, 0, 0 );
    }
  }

  return ans;
}

/**
 * Plotting of RooAbsPdf objects by generating the TGraphs using a RooPlot object
 * and plotting the graphs onto the TPad. The supported options are:
 *
 * - *Any* RooCmdArg supported by the RooAbsData::plotOn function. These will
 *   take precedence to custom defined options, with the excpetion of the
 *   RooFit::DrawOptions method.
 *
 * - PlotType: Specifying the representation for the data object on the Pad. Any
 *   of the types valid for the PlotGraph functions would be valid. If not
 *   specified, the simplefunc or fittedfunc methods will be used, depending on
 *   if the RooFit::VisualizeError is used.
 *
 * - EntryText: String to add in the legend entry. If this options is not
 *   present, then this object will NOT appear in the generated legend. Notice
 *   that the attributes to display in the legend would be generated from the
 *   PlotType used.
 *
 * - TrackY: Whether or not the y-axis range should be adjusted according to the
 *   newly added graph. By default, RooPdf objects would not be used to adjust
 *   the axis range.
 *
 * Additional automation of options include additional generation for
 * RooFit::VisualizeError. The stock RooPlot generates are contour line for the
 * uncertainty range rather than a line with error, making styling of a PDF with
 * uncertainty rather tedious. This functions takes the generated TGraphs by the
 * RooPlot and recalculated a TGraph with uncertainty. The newly calculated graph
 * will be placed under the ownership of the Pad.
 */
TGraph&
Pad1D::PlotPdf( RooAbsPdf& pdf, const std::vector<RooCmdArg>& arglist )
{
  const RooArgContainer args( arglist,
      {// Defining default arguments
        TrackY( TrackY::max ),// Default track y: only top
        RooArgContainer::CheckList( arglist, VisualizeError::CmdName ) ?
        PlotType( fittedfunc ) : PlotType( simplefunc )// default plot style.
      }
                              );
  TGraph& ans = MakePdfGraph( pdf, arglist );

  if( !args.Has( "Invisible" ) ){// Allow for invisible stuff to exist
    if( _prevrangetype == rangetype::aut ){
      // Forcing the range type to be histogram like.
      _prevrangetype = rangetype::hist;
    }
    PlotGraph( ans, args );
  }

  return ans;
}

/**
 * @brief Helper function for generating a TGraph representation of a RooAbsPdf
 * plot
 *
 * In the case the RooFit::Visualization is used, this function generates two
 * graphs via the in-built RooFit methods , one for the central plot, and the
 * other the contor of the uncertainty. This function then recalculates the two
 * graphs into a single TGraph with uncertainty, which is makes manipulating the
 * the plot results of a single RooAbsPdf instance more intuitive (changing the
 * plot style of a single object, rather than two).
 */
TGraph&
Pad1D::MakePdfGraph( RooAbsPdf& pdf, const RooArgContainer& args )
{
  // Augmenting the Visualize Error command. We have to manually insert a
  // RooArgSet containing all the floating variables directly in the
  // VisualizeError instance. Otherwise the generation of the error band might
  // fail.
  auto CorrectVisError
    = [this]( const RooAbsPdf& pdf,
              const usr::RooArgContainer& args ) -> RooCmdArg
      {
        if( args.Has<VisualizeError>() ){
          auto& cmd = args.Get<VisualizeError>();
          if( cmd.has_set() ){
            // If Visualized Parameter set is already specified, then simply
            // return the original parameter
            return cmd;
          } else {
            RooCmdArg ans( cmd );
            const auto& fit = cmd.GetRooFitResult();
            // Memory leak??
            RooArgSet* cloneParams = pdf.getObservables( fit.floatParsFinal() );
            ans.setSet( 0, *cloneParams );
            this->ClaimObject( cloneParams );
            return ans;
          }
        } else {
          return RooCmdArg::none();
        }
      };

  const RooCmdArg viscmd = CorrectVisError( pdf, args );
  RooLinkedList roolist  = args.MakeRooList( {VisualizeError::CmdName} );

  if( !args.Has<VisualizeError>() ){
    // Nothing special needs to be done for simple plotting.
    return GenGraph( pdf, roolist );
  } else {
    TGraph& centralgraph = GenGraph( pdf, roolist );
    roolist.Add( viscmd.Clone() );
    TGraph& errorgraph = GenGraph( pdf, roolist );

    // Map for storing uncertainty
    std::map<double, std::pair<double, double> > fiterr;

    for( int i = 0; i < errorgraph.GetN(); ++i ){
      const double x = errorgraph.GetX()[i];
      const double y = errorgraph.GetY()[i];
      if( !fiterr.count( x ) ){
        fiterr[x] = std::make_pair( 1e10, -1e10 );
      }
      fiterr.at( x ).first  = std::min( fiterr.at( x ).first, y  );
      fiterr.at( x ).second = std::max( fiterr.at( x ).second, y );
    }

    // The new TGraphAsymmErrors to store the results.
    TGraphAsymmErrors& ans = MakeObj<TGraphAsymmErrors>( fiterr.size() );
    unsigned i             = 0;// index for looping

    for( const auto& fiterrval : fiterr ){
      const double x   = fiterrval.first;
      const double y   = centralgraph.Eval( x );
      const double ylo = fiterrval.second.first;
      const double yhi = fiterrval.second.second;
      ans.SetPoint( i, x, y );
      ans.SetPointError( i, 0, 0, y-ylo, yhi-y );
      ++i;
    }

    return ans;
  }
}

// ------------------------------------------------------------------------------
// Private Helper functions
// ------------------------------------------------------------------------------

/**
 * @brief Helper function for generating RooAbsPdf plots onto a Pad.
 *
 * Throws and invalid_argument exception when the plotting call fails.
 */
TGraph&
Pad1D::GenGraph( RooAbsPdf& pdf, RooLinkedList& arglist )
{
  // Suppressing plotting messages
  RooMsgService::instance().setGlobalKillBelow( RooFit::WARNING );

  RooPlot* test = pdf.plotOn( &_frame, arglist );
  if( !test ){
    throw std::invalid_argument(
      "Bad argument list or object, plotting failed" );
  }

  auto& graph = _frame.LastPlot<TGraph>();
  graph.SetName( ( pdf.GetName() + RandomString( 6 ) ).c_str() );

  // Since this is a PDF and there is no weighting issue, we
  // enforce the fact that the graph should be positive definite.
  for( int i = 0; i < graph.GetN(); ++i ){
    const double x = graph.GetX()[i];
    const double y = graph.GetY()[i];
    if( y <= 0.0 ){
      graph.SetPoint( i, x, 1e-50 );
    }
  }

  return graph;
}

/**
 * @brief Helper function for generating RooAbsData plots onto a Pad.
 *
 * Throws an invalid_argument exception when the plotting call fails.
 */
TGraphAsymmErrors&
Pad1D::GenGraph( RooAbsData& data, RooLinkedList& arglist )
{
  // Suppressing plotting messages
  RooMsgService::instance().setGlobalKillBelow( RooFit::WARNING );

  // Generating plotting information
  RooPlot* test = data.plotOn( &_frame, arglist );
  if( !test ){
    throw std::invalid_argument(
      "Bad argument list or object, plotting failed" );
  }
  return _frame.LastPlot<TGraphAsymmErrors>();
}

/**
 * @brief Changing the stored _datamin, and _datamax variable according to object
 *
 * Moving to a private helper function to reduce verbosity in main implementation
 * function
 */
void
Pad1D::TrackObjectY( const TObject& obj, const int tracky )
{
  // Perfroming the axis range setting
  if( tracky == TrackY::min || tracky == TrackY::both ){
    if( obj.InheritsFrom( TH1D::Class() ) ){
      _datamin = std::min( _datamin,
        GetYmin( &dynamic_cast<const TH1D&>( obj ) ) );
    }
    if( obj.InheritsFrom( TGraph::Class() ) ){
      _datamin = std::min( _datamin,
        GetYmin( &dynamic_cast<const TGraph&>( obj ) ) );
    }
    if( obj.InheritsFrom( THStack::Class() ) ){
      _datamin = std::min( _datamin,
        GetYmin( &dynamic_cast<const THStack&>( obj ) ) );
    }
  }
  if( tracky == TrackY::max || tracky == TrackY::both ){
    if( obj.InheritsFrom( TH1D::Class() ) ){
      _datamax = std::max( _datamax,
        GetYmax( &dynamic_cast<const TH1D&>( obj ) ) );
    }
    if( obj.InheritsFrom( TGraph::Class() ) ){
      _datamax = std::max( _datamax,
        GetYmax( &dynamic_cast<const TGraph&>( obj ) ) );
    }
    if( obj.InheritsFrom( THStack::Class() ) ){
      _datamax = std::max( _datamax,
        GetYmax( &dynamic_cast<const THStack&>( obj ) ) );
    }
  }
  AutoSetYRange();
}

}/* plt */

}/* usr  */
