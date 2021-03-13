/**
 * @file    PlotCommon_Cmd.cc
 * @brief   Implementation of RooArgCmd related custom classes.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/RootUtils/RooArgContainer.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/Common/RootUtils/RooArgContainer.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include "RooFitResult.h"
#include "TFitResult.h"

#include <algorithm>
#include <set>
#include <string>

namespace usr {

namespace plt {

/**
 * @{
 * @addtogroup PlotUtilsArgument
 */
// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( PlotUnder );
/**
 * @{
 * @brief Request that the plotted object be placed underneath the target object
 **/
PlotUnder::PlotUnder( const TObject& obj ) :
  RooCmdArg( CmdName.c_str(),
             0, 0,// int
             0, 0,// double
             0, 0,// char
             &obj ){}

PlotUnder::PlotUnder( const TObject* obj ) :
  RooCmdArg( CmdName.c_str(),
             0, 0,// int
             0, 0,// double
             0, 0,// char
             obj ){}
/** @} */

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( TrackY );

/**
 * @brief define track y via enum
 */
TrackY::TrackY( const int i ) :
  RooCmdArg( CmdName.c_str(), i ){}

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( EntryText );

/**
 * @brief define entry text string.
 */
EntryText::EntryText( const std::string& str ) :
  RooCmdArg( CmdName.c_str(), 0, 0, 0, 0, str.c_str() ){}

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( PlotType );

/**
 * @brief defining plot type via enum
 */
PlotType::PlotType( const int i ) :
  RooCmdArg( CmdName.c_str(), i ){}

/**
 * @brief defining plot type via ROOT draw strings
 */
PlotType::PlotType( const std::string& drawopt ) :
  RooCmdArg( CmdName.c_str(), plottype_dummy, 0, 0, 0, drawopt.c_str() ){}

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( Plot2DF );

/**
 * @brief defining plot type via enum
 */
Plot2DF::Plot2DF( const int i ) :
  RooCmdArg( CmdName.c_str(), i ){}

/**
 * @brief defining plot type via ROOT draw strings
 */
Plot2DF::Plot2DF( const std::string& drawopt ) :
  RooCmdArg( CmdName.c_str(), plot2df_dummy, 0, 0, 0, drawopt.c_str() ){}

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( TextColor );

/**
 * @class TextColor
 * @brief Defining text color via the color index + alpha value convention.
 */
TextColor::TextColor( const int color, const float alpha ) :
  RooCmdArg( CmdName.c_str(), color, 0, alpha ){}

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( TextSize );

/**
 * @class TextSize
 * @brief Defining text size using the standard double notation (>1 is in
 * pixels), (<1 is in ratio to the pad width)
 */
TextSize::TextSize( const double size ) :
  RooCmdArg( CmdName.c_str(), 0, 0, size ){}

// ---------------------------------------------------------------------------//

USRUTILS_COMMON_REGISTERCMD( LineColor )
/**
 * @class LineColor
 * @brief Defining Line color via the color index + alpha value convention.
 */
LineColor::LineColor( const int color, const float alpha ) :
  RooCmdArg( CmdName.c_str(), color, 0, alpha ){}

USRUTILS_COMMON_REGISTERCMD( LineStyle );
/**
 * @class LineStyle
 * @brief Defining Line style via the style index (see namespace usr::plt::sty)
 */
LineStyle::LineStyle( const short style ) :
  RooCmdArg( CmdName.c_str(), style ){}

USRUTILS_COMMON_REGISTERCMD( LineWidth );
/**
 * @class LineWidth
 * @brief Defining line width via number of pixels
 */
LineWidth::LineWidth( const short width ) :
  RooCmdArg( CmdName.c_str(), width ){}

// ---------------------------------------------------------------------------//
USRUTILS_COMMON_REGISTERCMD( FillColor );
/**
 * @class FillColor
 * @brief Defining Fill color via color index and alpha value
 */
FillColor::FillColor( const int color, const float alpha ) :
  RooCmdArg( CmdName.c_str(), color, 0, alpha ){}

USRUTILS_COMMON_REGISTERCMD( FillStyle );
/**
 * @class LineStyle
 * @brief Defining fill style via the style index (see namespace usr::plt::sty)
 */
FillStyle::FillStyle( const short style ) :
  RooCmdArg( RooFit::FillStyle( style ) ){}

// ---------------------------------------------------------------------------//
USRUTILS_COMMON_REGISTERCMD( MarkerColor );
/**
 * @class MarkerColor
 * @brief Defining Marker color via color index and alpha value
 */
MarkerColor::MarkerColor( const int color, const float alpha ) :
  RooCmdArg( CmdName.c_str(), color, 0, alpha ){}

USRUTILS_COMMON_REGISTERCMD( MarkerStyle );
/**
 * @class MarkerStyle
 * @brief Defining marker style via the style index (see namespace usr::plt::sty)
 */
MarkerStyle::MarkerStyle( const short style ) :
  RooCmdArg( RooFit::MarkerStyle( style ) ){}

USRUTILS_COMMON_REGISTERCMD( MarkerSize );
/**
 * @class MarkerSize
 * @brief Defining MarkerSize  via number of pixels/8.
 */
MarkerSize::MarkerSize( const float size ) :
  RooCmdArg( RooFit::MarkerSize( size ) ){}

// ---------------------------------------------------------------------------//
USRUTILS_COMMON_REGISTERCMD( ScaleY );
ScaleY::ScaleY( const double s ) :
  RooCmdArg( CmdName.c_str(), 0, 0, s ){}


// ---------------------------------------------------------------------------//
USRUTILS_COMMON_REGISTERCMD( VisualizeError );
/**
 * @brief Switch for presenting fit uncertainty in fit.
 */
VisualizeError::VisualizeError( const TFitResultPtr& fit, const double z ) :
  RooCmdArg( CmdName.c_str(),
             0, 0,// int
             z, 0,// double
             0, 0,// c_string
             dynamic_cast<const TObject*>( &( *fit ) ) ){}

VisualizeError::VisualizeError( const RooFitResult& fit,
                                const double        z,
                                const bool          linearmethod ) :
  RooCmdArg( RooFit::VisualizeError( fit, z, linearmethod ) ){}

VisualizeError::VisualizeError( const RooFitResult& fit,
                                const RooArgSet&    param,
                                const double        z,
                                const bool          linearmethod ) :
  RooCmdArg( RooFit::VisualizeError( fit, param, z, linearmethod ) ){}

const RooFitResult&
VisualizeError::GetRooFitResult() const
{
  return *dynamic_cast<const RooFitResult*>( this->getObject( 0 ) );
}

const TFitResult&
VisualizeError::GetTFitResult() const
{
  return *dynamic_cast<const TFitResult*>( this->getObject( 0 ) );
}

const RooArgSet&
VisualizeError::set() const
{
  return *getSet( 0 );
};

bool
VisualizeError::has_set() const
{
  return getSet( 0 );
};

// ---------------------------------------------------------------------------//
USRUTILS_COMMON_REGISTERCMD( ExtrapolateInRatio );

ExtrapolateInRatio::ExtrapolateInRatio( const int flag ) :
  RooCmdArg( CmdName.c_str(), flag ){}


/** @} */

}

}
