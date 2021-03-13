/**
 * @file    Constants.cc
 * @brief   Defining the constants and extended constants
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#ifdef CMSSW_GIT_HASH
#include "UserUtils/PlotUtils/interface/Constants.hpp"
#else
#include "UserUtils/PlotUtils/Constants.hpp"
#endif

#include <vector>

namespace usr  {

namespace plt  {

/**
 * Plotting length in ROOT is defined using pixels (usually), with the
 * conversion of pixels into printable formats using the old standard of 72 DPI.
 *
 * Here we provide functions for convert from printing units (inches, mm... etc),
 * to the standard ROOT length units, as well as commonly used dimension in
 * CMS documents.
 */
namespace len {

extern const unsigned ROOT_DPI = 72;
extern const float INCH_TO_CM  = 2.54;

/**
 * @{
 * @brief Converting conventional units into ROOT standard length units.
 */
extern length_t
inch( const float x ){ return std::round( ROOT_DPI*x ); }

extern length_t
cm( const float x ){ return inch( x/INCH_TO_CM ); }

extern length_t
mm( const float x ){ return cm( x/10 ); }
/** @} */

extern unsigned
a4paperwidth(){ return mm( 210 ); }

extern unsigned
a4paperheight(){ return mm( 297 ); }

/**
 * @{
 * @brief Document dimensions taken from cms-tdr.cls
 */
extern unsigned
a4textwidth_default(){ return mm( 160 ); }

extern unsigned
a4textheight_default(){ return mm( 235 ); }
/** @} */

}/* len */


/**
 * @class usr::plt::font
 * Since font setting are more commonly copied around, implemenation of font
 * conversion functions would be a class to enable storage of variables.
 *
 * Fonts in ROOT (when using precision greater than 3) are defined in pixels,
 * which, when used together with the ROOT DPI of 72, makes it the common
 * unit of pts.
 */

/** @brief default font setting for CMS documents */
float font::default_basesize = 11;

/** @brief Recommended plot font for CMS documents */
font::face font::default_face = font::helvetica;

/** @brief default font line spacing multiplier */
float font::default_lineheight = 1.2;

/**
 * @brief Defining a font setting.
 * @param basesize   The base size of the font
 * @param face       The type-face of the font to use.
 * @param lineheight The lineheight **multiplier** to use.
 */
font::font( float basesize,
            face  face,
            float lineheight ) :
  _basesize( basesize ),
  _face( face ),
  _lineheight( em( lineheight ) )
{
}

/** @brief returning base size */
float
font::size() const { return _basesize; }

/** @brief returning line height (absolute, not relative to basesize) */
float
font::lineheight() const { return _lineheight; }

/** @brief returning the font face as standard ROOT font type */
short
font::fontface() const
{
  return ( _face * 10 )/10 + 3;
}

/** @brief returning scaled font size */
float
font::em( const float x ) const
{
  return x * _basesize;
}

/*-----------------------------------------------------------------------------
 *  Font-size conversion helper variables.
   --------------------------------------------------------------------------*/
static const float scaletable[7][4] = {
  // 10      11          12             others
  {
    5/10,     6/11.,      6/12.,         0.5// tiny
  },{
    8/10,     9/11.,      10/12.,        0.75// footnote
  },{
    9/10,     10/11.,     10.95/12.,     0.9// small
  },{
    12/10,    12/11.,     14.4/12.,      1.2// large
  },{
    14.4/10,  14.4/11.,   17.28/12.,     14.4// Large
  },{
    17.28/10, 17.28/11.,  20.74/12.,     1.72// LARGE
  },{
    20.7/10,  20.7/11.,   24.8/12.,      2   }// huge
};

enum e_size : short
{
  e_tiny,
  e_footnote,
  e_small,
  e_large,
  e_llarge,
  e_lllarge,
  e_huge
};

static unsigned
getcol( const float x )
{
  return x == 10 ? 0 :
         x == 11 ? 1 :
         x == 12 ? 2 :
         3;
}

/** @{
 * @brief conversion to LateX [in-build font sizes](https://en.wikibooks.org/wiki/LaTeX/Fonts#Built-in_sizes)
 */
float
font::tiny() const
{ return _basesize*scaletable[e_tiny][getcol( _basesize )]; }

float
font::footnote() const
{ return _basesize*scaletable[e_footnote][getcol( _basesize )]; }

float
font::small() const
{ return _basesize*scaletable[e_small][getcol( _basesize )]; }

float
font::large() const
{ return _basesize*scaletable[e_large][getcol( _basesize )]; }

float
font::Large() const
{ return _basesize*scaletable[e_llarge][getcol( _basesize )]; }

float
font::LARGE() const
{ return _basesize*scaletable[e_lllarge][getcol( _basesize )]; }

float
font::huge() const
{ return _basesize*scaletable[e_huge][getcol( _basesize )]; }

/** @} */


/**
 * Notice that these strings are to be used in the ROOT plotting functions,
 * so it will be using the ROOT style '#' for latex commands rather than
 * backslashes
 */
namespace unit  {
const std::string GeV   = "GeV";
const std::string GeVc  = "GeV/#it{c}";
const std::string GeVcc = "GeV/#it{c}^{2}";
}/* unit  */


namespace cap  {
const std::string prelim = "Preliminary";
const std::string sim    = "Simulation";
}/* cap  */

/**
 * @details  Probably one of the worst offenders of magic-number usage in @ROOT:
 * that of the styling variables. Unless you have memorised the contents of the
 * documentation, you would probably have no idea what someone else code is doing
 * to stylise the plots.
 */
namespace sty {

/**
 * @{
 * @brief Predefined fill styles for @ROOT{TAttrFill}.
 * */
const short fillnone      = 0;
const short fillsolid     = 1001;
const short filldotdense  = 3001;
const short filldot       = 3002;
const short filldotsparse = 3003;
/** @} */

/**
 * @brief  A human readable interface to hash-fill styles in @ROOT{TAttFill}.
 *
 * To be used with the distmm(), angle1(), and angle2() functions to make what
 * the hash is suppose to look like look more natural. The call to a function
 * would SetFillStyle Function would then be something like:
 * ```
 * myobject.SetFillStyle( fillhash(distmm(0.9), angle1(90), angle2(100)));
 * ```
 *
 * @param  distance distance between hash lines
 * @param  angle1   angle of first hash lines
 * @param  angle2   angle of second hash lines
 * @return          the value corresponding to the angle code.
 */
const short
fillhash( unsigned short distance,
          unsigned short angle1,
          unsigned short angle2 )
{
  distance = std::max( distance, (unsigned short)1 );
  distance = std::min( distance, (unsigned short)9 );
  angle1   = std::min( angle1,  (unsigned short)9 );
  angle2   = std::min( angle2,  (unsigned short)9 );
  return 3000 + 100*distance+10*angle1+1*angle2;
}

/**
 * @brief Converting actual hash distances in millimeters into the unit
 * used by the distance string.
 *
 * @ROOT{TAttFill} only defines 0.5mm converts to 1 , and 6mm converts to 9.
 * There we assume a linear relation, rounding to the closest value.
 */
const short
distmm( const float x )
{
  return x < usr::plt::len::mm(0.84375) ? 1 :
         x < usr::plt::len::mm(1.53125) ? 2 :
         x < usr::plt::len::mm(2.21875) ? 3 :
         x < usr::plt::len::mm(2.90625) ? 4 :
         x < usr::plt::len::mm(3.59375) ? 5 :
         x < usr::plt::len::mm(4.28125) ? 6 :
         x < usr::plt::len::mm(4.96875) ? 7 :
         x < usr::plt::len::mm(5.65625) ? 8 :
         9;
}

/**
 * @brief Converting the first angle of the hash in degrees into defined digits.
 *
 * Angles are defined in 10 degrees increments (with the exception of the 45
 * degree mark), rounding to closest defined value.
 */
const short
angle1( const float x )
{
  return x < 0    ? anglenone :
         x < 5    ? 0 :
         x < 15   ? 1 :
         x < 25   ? 2 :
         x < 37.5 ? 3 :
         x < 52.5 ? 4 :
         x < 65   ? 6 :
         x < 75   ? 7 :
         x < 85   ? 8 :
         x <= 90  ? 9 :
         anglenone;
}

/**
 * @brief Converting the first angle of the hash in degrees into defined digits.
 *
 * Angles are defined in 10 degrees increments (with the exception of the 45
 * degree mark), rounding to closest defined value.
 */
const short
angle2( const float x )
{
  return x > 180   ? anglenone :
         x > 175   ? 0 :
         x > 165   ? 1 :
         x > 155   ? 2 :
         x > 142.5 ? 3 :
         x > 127.5 ? 4 :
         x > 115   ? 6 :
         x > 105   ? 7 :
         x > 95    ? 8 :
         x >= 90   ? 9 :
         anglenone;
}

/** @brief Flag to use if you want to diplay a certain angle display */
const short anglenone = 5;

/**
 * @{
 * @brief common settings used by @ROOT{TAttLine}.
 * */
const short linsolid     = 1;
const short lindotted    = 2;
const short lindensedot  = 3;
const short lindashed    = 9;
const short linshortdash = 7;
/** @} */

/**
 * @{
 * @brief common settings used by @ROOT{TAttMarker}.
 * */
extern const short mkrcircle           = 20;
extern const short mkrsquare           = 21;
extern const short mkrtriangleup       = 22;
extern const short mkrtriangledown     = 23;
extern const short mkrdiamond          = 33;
extern const short mkrcross            = 34;
extern const short mkrstar             = 29;
extern const short mkropencircle       = 24;
extern const short mkropensquare       = 25;
extern const short mkropentriangleup   = 26;
extern const short mkropentriangledown = 32;
extern const short mkropendiamond      = 27;
extern const short mkropencross        = 28;
extern const short mkropenstar         = 30;
/** @} */

}/* sty */

}/* plt  */

}/* usr  */
