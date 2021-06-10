#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#else
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#endif

namespace usr {

namespace plt {

/**
 * @class RangeByVar
 * Constructor containter for constructing a Pad1D instance by using RooRealVar
 * as the x axis. By passing this class into the Pad1D constructor, the Pad1D
 * would know to generate a corresponding RooPlot instance for RooFit plot
 * object generation.
 */
const double RangeByVar::find_default
  = std::numeric_limits<double>::lowest();

/**
 * @brief Define a range using the @ROOT{RooRealVar} reference.
 *
 * Could defined x range and number of bins to use for the x axis.
 */
RangeByVar::RangeByVar(
  const RooRealVar& _var,
  const double      _xmin,
  const double      _xmax,
  const int         _nbin
  ) : var( &_var )
{
  xmin = _xmin == find_default ? var->getMin() :
         _xmin > _xmax         ? var->getMin() :
         _xmin;
  xmax = _xmax == find_default ? var->getMax() :
         _xmin < _xmax         ? var->getMax() :
         _xmax;
  nbin = _nbin < 0 ? var->getBins() : _nbin;
}

/**
 * @brief Define a range using the @ROOT{RooRealVar} pointer
 *
 * Could defined x range and number of bins to use for the x axis.
 */
RangeByVar::RangeByVar(
  const RooRealVar* _var,
  const double      _xmin,
  const double      _xmax,
  const int         _nbin
  ) : var( _var )
{
  xmin = _xmin == find_default ? var->getMin() :
         _xmin > _xmax         ? var->getMin() :
         _xmin;
  xmax = _xmax == find_default ? var->getMax() :
         _xmin > _xmax         ? var->getMax() :
         _xmax;
  nbin = _nbin < 0 ? var->getBins() : _nbin;
}


/**
 * @class RooFrame
 * @ingroup PlotUtilsCore
 * @details
 * An simple wrapper object designed to expose additional features of the
 * RooPlot class (the internal TH1 used for axis plotting) and a feature for
 * allow the creation objects directly under the ownership of a RooPlot instance
 */

/**
 * @brief Construct a frame with a variable.
 * @brief Name of the object would be assigned a random string to avoid
 * collision
 */
RooFrame::RooFrame( const RangeByVar& range ) :
  RooPlot( RandomString( 12 ).c_str(), "",
           *range.var, range.xmin, range.xmax, range.nbin )
{

}

/**
 * @brief Construct a frame with just the x axis range.
 */
RooFrame::RooFrame( const double min, const double max ) :
  RooPlot( min, max )
{
  RooPlot::SetName( RandomString( 12 ).c_str() );
  RooPlot::SetTitle( "" );
}

/**
 * @brief Construct a frame with just the x axis range.
 */
RooFrame::RooFrame( const double xmin, const double xmax,
                    const double ymin, const double ymax ) :
  RooPlot( xmin, xmax, ymin, ymax )
{
  RooPlot::SetName( RandomString( 12 ).c_str() );
  RooPlot::SetTitle( "" );
}

/**
 * @brief Must allow for null constructor for pad that doesn't use a frame.
 */
RooFrame::RooFrame() :
  RooPlot()
{// Do nothing
}

RooFrame::~RooFrame(){}

}

}