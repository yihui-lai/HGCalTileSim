/**
 * @file    Constants.hpp
 * @brief   Unit conversion helper functions
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifndef USERUTILS_PLOTUTILS_UNITS_HPP
#define USERUTILS_PLOTUTILS_UNITS_HPP

#include <cmath>
#include <string>

namespace usr  {

namespace plt  {

/**
 * @defgroup PlotUtilsConstants PlotUtilsConstants
 * @ingroup PlotUtils
 * @brief ROOT plotting's magic number conversion functions and constants.
 * @{
 */

/** @brief Abstact length unit used in the plotting routines in root */
typedef unsigned length_t;

/**
 * @brief length conversion functions and constants.
 */
namespace len {

/**
 * @addtogroup PlotUtilsConstants
 * @{
 */
extern const unsigned ROOT_DPI;

extern length_t inch( const float x );
extern length_t cm( const float x );
extern length_t mm( const float x );

extern length_t a4paperwidth();
extern length_t a4paperheight();
extern length_t a4textwidth_default();
extern length_t a4textheight_default();
/** @} */

}/* len */

/**
 * @brief Container for a font setting and size conversion.
 * @ingroup PlotUtilsConstants
 */
class font
{
public:// static and user enums.

  /** @brief enum interface for defined fonts. */
  enum face : short
  {
    helvetica = 40,
    courier   = 80,
    times     = 130
  };

  /** @brief enum interface for text alignment in @ROOT{TLatex} */
  enum align : short
  {
    left          = 10,
    hcenter       = 20,
    right         = 30,
    top           = 3,
    vcenter       = 2,
    bottom        = 1,
    bottom_left   = 11,
    bottom_center = 21,
    bottom_right  = 31,
    center_left   = 12,
    center        = 22,
    center_right  = 32,
    top_left      = 13,
    top_center    = 23,
    top_right     = 33
  };

  static float default_basesize;
  static font::face default_face;
  static float default_lineheight;

public:
  font( float basesize   = default_basesize,
        face  face       = default_face,
        float lineheight = default_lineheight );

  float size() const;
  float lineheight() const;
  short fontface() const;
  float em( const float mult ) const;

  float tiny() const;
  float footnote() const;
  float small() const;
  float large() const;
  float Large() const;
  float LARGE() const;
  float huge() const;

protected:
  float _basesize;
  font::face _face;
  float _lineheight;
};

using align = font::align;// Align can be used in more places than just text.

/** @brief Strings for a unified presentation of units in plots */
namespace unit {
extern const std::string GeV;
extern const std::string GeVc;
extern const std::string GeVcc;
}/* unit */

/** @brief String for unifed presentation of common captions text in plots */
namespace cap  {
extern const std::string prelim;
extern const std::string sim;
}/* cap  */

/** @brief root style conversion functions  */
namespace sty  {

/**
 * @addtogroup PlotUtilsConstants
 * @{
 */

// Marker styles
inline float
mkrpixels( const int x ){ return x / 8; }

extern const short mkrcircle;
extern const short mkrsquare;
extern const short mkrtriangleup;
extern const short mkrtriangledown;
extern const short mkrdiamond;
extern const short mkrcross;
extern const short mkrstar;
extern const short mkropencircle;
extern const short mkropensquare;
extern const short mkropentriangleup;
extern const short mkropentriangledown;
extern const short mkropendiamond;
extern const short mkropencross;
extern const short mkropenstar;


extern const short fillnone;
extern const short fillsolid;
extern const short filldotdense;
extern const short filldot;
extern const short filldotsparse;

const short fillhash( unsigned short distance,
                      unsigned short angle1,
                      unsigned short angle2 );
extern const short distmm( const float );
extern const short anglenone;
extern const short angle1( const float );
extern const short angle2( const float );

// Line style
extern const short linsolid;
extern const short lindotted;
extern const short lindensedot;
extern const short lindashed;
extern const short linshortdash;

/** @} */

}/* sty  */

/**
 * @brief Defining a more comprehensive list of colors
 */
namespace col
{

/**
 * @addtogroup PlotUtilsConstants
 * @{
 */
extern const int maroon;
extern const int darkred;
extern const int brown;
extern const int firebrick;
extern const int crimson;
extern const int red;
extern const int tomato;
extern const int coral;
extern const int indianred;
extern const int lightcoral;
extern const int darksalmon;
extern const int salmon;
extern const int lightsalmon;
extern const int orangered;
extern const int darkorange;
extern const int orange;
extern const int gold;
extern const int darkgoldenrod;
extern const int goldenrod;
extern const int palegoldenrod;
extern const int darkkhaki;
extern const int khaki;
extern const int olive;
extern const int yellow;
extern const int yellowgreen;
extern const int darkolivegreen;
extern const int olivedrab;
extern const int lawngreen;
extern const int chartreuse;
extern const int greenyellow;
extern const int darkgreen;
extern const int green;
extern const int forestgreen;
extern const int lime;
extern const int limegreen;
extern const int lightgreen;
extern const int palegreen;
extern const int darkseagreen;
extern const int mediumspringgreen;
extern const int springgreen;
extern const int seagreen;
extern const int mediumaquamarine;
extern const int mediumseagreen;
extern const int lightseagreen;
extern const int darkslategray;
extern const int teal;
extern const int darkcyan;
extern const int aqua;
extern const int cyan;
extern const int lightcyan;
extern const int darkturquoise;
extern const int turquoise;
extern const int mediumturquoise;
extern const int paleturquoise;
extern const int aquamarine;
extern const int powderblue;
extern const int cadetblue;
extern const int steelblue;
extern const int cornflowerblue;
extern const int deepskyblue;
extern const int dodgerblue;
extern const int lightblue;
extern const int skyblue;
extern const int lightskyblue;
extern const int midnightblue;
extern const int navy;
extern const int darkblue;
extern const int mediumblue;
extern const int blue;
extern const int royalblue;
extern const int blueviolet;
extern const int indigo;
extern const int darkslateblue;
extern const int slateblue;
extern const int mediumslateblue;
extern const int mediumpurple;
extern const int darkmagenta;
extern const int darkviolet;
extern const int darkorchid;
extern const int mediumorchid;
extern const int purple;
extern const int thistle;
extern const int plum;
extern const int violet;
extern const int magenta;
extern const int fuchsia;
extern const int orchid;
extern const int mediumvioletred;
extern const int palevioletred;
extern const int deeppink;
extern const int hotpink;
extern const int lightpink;
extern const int pink;
extern const int antiquewhite;
extern const int beige;
extern const int bisque;
extern const int blanchedalmond;
extern const int wheat;
extern const int cornsilk;
extern const int lemonchiffon;
extern const int lightgoldenrodyellow;
extern const int lightyellow;
extern const int saddlebrown;
extern const int sienna;
extern const int chocolate;
extern const int peru;
extern const int sandybrown;
extern const int burlywood;
extern const int tan;
extern const int rosybrown;
extern const int moccasin;
extern const int navajowhite;
extern const int peachpuff;
extern const int mistyrose;
extern const int lavenderblush;
extern const int linen;
extern const int oldlace;
extern const int papayawhip;
extern const int seashell;
extern const int mintcream;
extern const int slategray;
extern const int lightslategray;
extern const int lightsteelblue;
extern const int lavender;
extern const int floralwhite;
extern const int aliceblue;
extern const int ghostwhite;
extern const int honeydew;
extern const int ivory;
extern const int azure;
extern const int snow;
extern const int black;
extern const int dimgrey;
extern const int grey;
extern const int darkgrey;
extern const int silver;
extern const int lightgrey;
extern const int gainsboro;
extern const int whitesmoke;
extern const int white;
extern const int fuchsia;
extern const int dimgray;
extern const int darkgray;
extern const int lightgray;
extern const int gray;
extern const int bzgreen;
extern const int bzyellow;
/** @} */

}// col

/** @} */

}/* plt  */

}/* usr  */

#endif/* end of include guard: USERUTILS_PLOTUTILS_CANVAS_HPP */
