#include "TColor.h"

namespace usr {

namespace plt {

/**
 * @details Most of the colors are defined in:
 * https://www.rapidtables.com/web/color/RGB_Color.html#color-table
 *
 * Which includes both British-English and American-English conventions of
 * gray/grey. A few more standard colors used for CMS plots have also been
 * defined.
 */
namespace col {

extern const int maroon               = TColor::GetColor( "#800000" );
extern const int darkred              = TColor::GetColor( "#8B0000" );
extern const int brown                = TColor::GetColor( "#A52A2A" );
extern const int firebrick            = TColor::GetColor( "#B22222" );
extern const int crimson              = TColor::GetColor( "#DC143C" );
extern const int red                  = TColor::GetColor( "#FF0000" );
extern const int tomato               = TColor::GetColor( "#FF6347" );
extern const int coral                = TColor::GetColor( "#FF7F50" );
extern const int indianred            = TColor::GetColor( "#CD5C5C" );
extern const int lightcoral           = TColor::GetColor( "#F08080" );
extern const int darksalmon           = TColor::GetColor( "#E9967A" );
extern const int salmon               = TColor::GetColor( "#FA8072" );
extern const int lightsalmon          = TColor::GetColor( "#FFA07A" );
extern const int orangered            = TColor::GetColor( "#FF4500" );
extern const int darkorange           = TColor::GetColor( "#FF8C00" );
extern const int orange               = TColor::GetColor( "#FFA500" );
extern const int gold                 = TColor::GetColor( "#FFD700" );
extern const int darkgoldenrod        = TColor::GetColor( "#B8860B" );
extern const int goldenrod            = TColor::GetColor( "#DAA520" );
extern const int palegoldenrod        = TColor::GetColor( "#EEE8AA" );
extern const int darkkhaki            = TColor::GetColor( "#BDB76B" );
extern const int khaki                = TColor::GetColor( "#F0E68C" );
extern const int olive                = TColor::GetColor( "#808000" );
extern const int yellow               = TColor::GetColor( "#FFFF00" );
extern const int yellowgreen          = TColor::GetColor( "#9ACD32" );
extern const int darkolivegreen       = TColor::GetColor( "#556B2F" );
extern const int olivedrab            = TColor::GetColor( "#6B8E23" );
extern const int lawngreen            = TColor::GetColor( "#7CFC00" );
extern const int chartreuse           = TColor::GetColor( "#7FFF00" );
extern const int greenyellow          = TColor::GetColor( "#ADFF2F" );
extern const int darkgreen            = TColor::GetColor( "#006400" );
extern const int green                = TColor::GetColor( "#008000" );
extern const int forestgreen          = TColor::GetColor( "#228B22" );
extern const int lime                 = TColor::GetColor( "#00FF00" );
extern const int limegreen            = TColor::GetColor( "#32CD32" );
extern const int lightgreen           = TColor::GetColor( "#90EE90" );
extern const int palegreen            = TColor::GetColor( "#98FB98" );
extern const int darkseagreen         = TColor::GetColor( "#8FBC8F" );
extern const int mediumspringgreen    = TColor::GetColor( "#00FA9A" );
extern const int springgreen          = TColor::GetColor( "#00FF7F" );
extern const int seagreen             = TColor::GetColor( "#2E8B57" );
extern const int mediumaquamarine     = TColor::GetColor( "#66CDAA" );
extern const int mediumseagreen       = TColor::GetColor( "#3CB371" );
extern const int lightseagreen        = TColor::GetColor( "#20B2AA" );
extern const int darkslategray        = TColor::GetColor( "#2F4F4F" );
extern const int teal                 = TColor::GetColor( "#008080" );
extern const int darkcyan             = TColor::GetColor( "#008B8B" );
extern const int aqua                 = TColor::GetColor( "#00FFFF" );
extern const int cyan                 = TColor::GetColor( "#00FFFF" );
extern const int lightcyan            = TColor::GetColor( "#E0FFFF" );
extern const int darkturquoise        = TColor::GetColor( "#00CED1" );
extern const int turquoise            = TColor::GetColor( "#40E0D0" );
extern const int mediumturquoise      = TColor::GetColor( "#48D1CC" );
extern const int paleturquoise        = TColor::GetColor( "#AFEEEE" );
extern const int aquamarine           = TColor::GetColor( "#7FFFD4" );
extern const int powderblue           = TColor::GetColor( "#B0E0E6" );
extern const int cadetblue            = TColor::GetColor( "#5F9EA0" );
extern const int steelblue            = TColor::GetColor( "#4682B4" );
extern const int cornflowerblue       = TColor::GetColor( "#6495ED" );
extern const int deepskyblue          = TColor::GetColor( "#00BFFF" );
extern const int dodgerblue           = TColor::GetColor( "#1E90FF" );
extern const int lightblue            = TColor::GetColor( "#ADD8E6" );
extern const int skyblue              = TColor::GetColor( "#87CEEB" );
extern const int lightskyblue         = TColor::GetColor( "#87CEFA" );
extern const int midnightblue         = TColor::GetColor( "#191970" );
extern const int navy                 = TColor::GetColor( "#000080" );
extern const int darkblue             = TColor::GetColor( "#00008B" );
extern const int mediumblue           = TColor::GetColor( "#0000CD" );
extern const int blue                 = TColor::GetColor( "#0000FF" );
extern const int royalblue            = TColor::GetColor( "#4169E1" );
extern const int blueviolet           = TColor::GetColor( "#8A2BE2" );
extern const int indigo               = TColor::GetColor( "#4B0082" );
extern const int darkslateblue        = TColor::GetColor( "#483D8B" );
extern const int slateblue            = TColor::GetColor( "#6A5ACD" );
extern const int mediumslateblue      = TColor::GetColor( "#7B68EE" );
extern const int mediumpurple         = TColor::GetColor( "#9370DB" );
extern const int darkmagenta          = TColor::GetColor( "#8B008B" );
extern const int darkviolet           = TColor::GetColor( "#9400D3" );
extern const int darkorchid           = TColor::GetColor( "#9932CC" );
extern const int mediumorchid         = TColor::GetColor( "#BA55D3" );
extern const int purple               = TColor::GetColor( "#800080" );
extern const int thistle              = TColor::GetColor( "#D8BFD8" );
extern const int plum                 = TColor::GetColor( "#DDA0DD" );
extern const int violet               = TColor::GetColor( "#EE82EE" );
extern const int magenta              = TColor::GetColor( "#FF00FF" );
extern const int fuchsia              = TColor::GetColor( "#FF00FF" );
extern const int orchid               = TColor::GetColor( "#DA70D6" );
extern const int mediumvioletred      = TColor::GetColor( "#C71585" );
extern const int palevioletred        = TColor::GetColor( "#DB7093" );
extern const int deeppink             = TColor::GetColor( "#FF1493" );
extern const int hotpink              = TColor::GetColor( "#FF69B4" );
extern const int lightpink            = TColor::GetColor( "#FFB6C1" );
extern const int pink                 = TColor::GetColor( "#FFC0CB" );
extern const int antiquewhite         = TColor::GetColor( "#FAEBD7" );
extern const int beige                = TColor::GetColor( "#F5F5DC" );
extern const int bisque               = TColor::GetColor( "#FFE4C4" );
extern const int blanchedalmond       = TColor::GetColor( "#FFEBCD" );
extern const int wheat                = TColor::GetColor( "#F5DEB3" );
extern const int cornsilk             = TColor::GetColor( "#FFF8DC" );
extern const int lemonchiffon         = TColor::GetColor( "#FFFACD" );
extern const int lightgoldenrodyellow = TColor::GetColor( "#FAFAD2" );
extern const int lightyellow          = TColor::GetColor( "#FFFFE0" );
extern const int saddlebrown          = TColor::GetColor( "#8B4513" );
extern const int sienna               = TColor::GetColor( "#A0522D" );
extern const int chocolate            = TColor::GetColor( "#D2691E" );
extern const int peru                 = TColor::GetColor( "#CD853F" );
extern const int sandybrown           = TColor::GetColor( "#F4A460" );
extern const int burlywood            = TColor::GetColor( "#DEB887" );
extern const int tan                  = TColor::GetColor( "#D2B48C" );
extern const int rosybrown            = TColor::GetColor( "#BC8F8F" );
extern const int moccasin             = TColor::GetColor( "#FFE4B5" );
extern const int navajowhite          = TColor::GetColor( "#FFDEAD" );
extern const int peachpuff            = TColor::GetColor( "#FFDAB9" );
extern const int mistyrose            = TColor::GetColor( "#FFE4E1" );
extern const int lavenderblush        = TColor::GetColor( "#FFF0F5" );
extern const int linen                = TColor::GetColor( "#FAF0E6" );
extern const int oldlace              = TColor::GetColor( "#FDF5E6" );
extern const int papayawhip           = TColor::GetColor( "#FFEFD5" );
extern const int seashell             = TColor::GetColor( "#FFF5EE" );
extern const int mintcream            = TColor::GetColor( "#F5FFFA" );
extern const int slategray            = TColor::GetColor( "#708090" );
extern const int lightslategray       = TColor::GetColor( "#778899" );
extern const int lightsteelblue       = TColor::GetColor( "#B0C4DE" );
extern const int lavender             = TColor::GetColor( "#E6E6FA" );
extern const int floralwhite          = TColor::GetColor( "#FFFAF0" );
extern const int aliceblue            = TColor::GetColor( "#F0F8FF" );
extern const int ghostwhite           = TColor::GetColor( "#F8F8FF" );
extern const int honeydew             = TColor::GetColor( "#F0FFF0" );
extern const int ivory                = TColor::GetColor( "#FFFFF0" );
extern const int azure                = TColor::GetColor( "#F0FFFF" );
extern const int snow                 = TColor::GetColor( "#FFFAFA" );
extern const int black                = TColor::GetColor( "#000000" );
extern const int dimgray              = TColor::GetColor( "#696969" );
extern const int dimgrey              = TColor::GetColor( "#696969" );
extern const int gray                 = TColor::GetColor( "#808080" );
extern const int grey                 = TColor::GetColor( "#808080" );
extern const int darkgray             = TColor::GetColor( "#A9A9A9" );
extern const int darkgrey             = TColor::GetColor( "#A9A9A9" );
extern const int silver               = TColor::GetColor( "#C0C0C0" );
extern const int lightgray            = TColor::GetColor( "#D3D3D3" );
extern const int lightgrey            = TColor::GetColor( "#D3D3D3" );
extern const int gainsboro            = TColor::GetColor( "#DCDCDC" );
extern const int whitesmoke           = TColor::GetColor( "#F5F5F5" );
extern const int white                = TColor::GetColor( "#FFFFFF" );

// Standard colors for Brazillean flag plot
extern const int bzgreen  = kGreen + 1;
extern const int bzyellow = kOrange;

}

}

}
