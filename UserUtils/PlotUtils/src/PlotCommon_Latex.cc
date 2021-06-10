/**
 * @file PlotCommon_Latex.cc
 * @author Yi-Mu "Enoch" Chen
 *
 * @brief Computation of Latex dimension estimates.
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/interface/PlotCommon.hpp"
#else
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/PlotUtils/PlotCommon.hpp"
#endif

#include <map>
#include <regex>

static const std::regex dec_regex( ".*?(#[a-zA-Z]+\\b\\s*\\{).*?" );
static const std::regex sym_regex( "#[a-zA-Z]+\\b" );
static const double subratio       = 2./3.;
static const double fraclineheight = 1./5.;
static const double int_height     = 2.0;

// Precalculated character width chart (relative to font size)
static const std::map<char, double> char_width = {
  {' ', 0.30120482}, {'!', 0.20080321}, {'$', 0.50200803}, {'%', 0.80321285},
  {'&', 0.60240964}, {'(', 0.30120482}, {')', 0.25100402}, {'*', 0.35140562},
  {'+', 0.50200803}, {',', 0.20080321}, {'-', 0.30120482}, {'.', 0.20080321},
  {'/', 0.35140562}, {'0', 0.50200803}, {'1', 0.35140562}, {'2', 0.50200803},
  {'3', 0.50200803}, {'4', 0.50200803}, {'5', 0.50200803}, {'6', 0.50200803},
  {'7', 0.50200803}, {'8', 0.50200803}, {'9', 0.50200803}, {';', 0.20080321},
  {'<', 0.50200803}, {'=', 0.50200803}, {'>', 0.50200803}, {'@', 0.90361446},
  {'A', 0.60240964}, {'B', 0.60240964}, {'C', 0.65261044}, {'D', 0.65261044},
  {'E', 0.60240964}, {'F', 0.55220884}, {'G', 0.65261044}, {'H', 0.60240964},
  {'I', 0.20080321}, {'J', 0.40160643}, {'K', 0.60240964}, {'L', 0.50200803},
  {'M', 0.70281124}, {'N', 0.60240964}, {'O', 0.70281124}, {'P', 0.55220884},
  {'Q', 0.70281124}, {'R', 0.65261044}, {'S', 0.60240964}, {'T', 0.55220884},
  {'U', 0.60240964}, {'V', 0.60240964}, {'W', 0.85341365}, {'X', 0.60240964},
  {'Y', 0.60240964}, {'Z', 0.55220884}, {'[', 0.25100402}, {']', 0.20080321},
  {'a', 0.50200803}, {'b', 0.50200803}, {'c', 0.45180723}, {'d', 0.45180723},
  {'e', 0.50200803}, {'f', 0.25100402}, {'g', 0.45180723}, {'h', 0.45180723},
  {'i', 0.15060241}, {'j', 0.20080321}, {'k', 0.45180723}, {'l', 0.15060241},
  {'m', 0.70281124}, {'n', 0.45180723}, {'o', 0.50200803}, {'p', 0.50200803},
  {'q', 0.45180723}, {'r', 0.30120482}, {'s', 0.45180723}, {'t', 0.25100402},
  {'u', 0.45180723}, {'v', 0.45180723}, {'w', 0.65261044}, {'x', 0.45180723},
  {'y', 0.45180723}, {'z', 0.45180723}
};

namespace usr {

namespace plt {

/**
 * @brief Better estimation of the actual width (in units 'em') of a Latex
 * string.
 *
 * @details This functions sums the width of most printable ASCII characters with
 * by a predefined table, as the Latex engine used in @(ROOT} has very minimal
 * kerning on the characeter sets. The table itself is generated using a ROOT
 * script, using a Helvetic font (so estimates for Times font might be off).
 * Special care is taken to remove decorators such as fractions, square-roots,
 * left-right braces, and sub/sup-scripts. Non-ASCII characters (greek symbols,
 * mathematical symbols... etc), are simply estimated to be a standard 'M' width.
 */
double
EstimateLatexWidth( const std::string& text )
{
  double ans = 0;
  std::string un_text( text );
  std::smatch match;
  size_t find = std::string::npos;

  // Dealing with fractions in text
  while( ( find = un_text.find( "#frac{" ) ) != std::string::npos ){
    const size_t openbrace1  = NextOpenBrace( un_text, find );
    const size_t closebrace1 = MatchBrace( un_text, openbrace1 );
    const size_t openbrace2  = NextOpenBrace( un_text, closebrace1 );
    const size_t closebrace2 = MatchBrace( un_text, openbrace2 );

    const std::string numstr( un_text.begin() + openbrace1+1,
                              un_text.begin() + closebrace1 );
    const std::string denstr( un_text.begin() + openbrace2+1,
                              un_text.begin() + closebrace2 );
    ans += std::max( EstimateLatexWidth( numstr ),
      EstimateLatexWidth( denstr ) );

    un_text.erase( un_text.begin() + find, un_text.begin() + closebrace2+1 );
  }

  // Dealing with square roots
  while( ( find = un_text.find( "#sqrt" ) ) != std::string::npos ){
    const size_t openbrace  = NextOpenBrace( un_text, find );
    const size_t closebrace = MatchBrace( un_text, openbrace );

    un_text.erase( un_text.begin() + closebrace );
    un_text.erase( un_text.begin() + find, un_text.begin() + openbrace + 1 );
    ans += 1;
  }

  // Dealing with super scripts and sub_scripts
  while( ( find = un_text.find( "^{" ) ) != std::string::npos ){
    const size_t openbrace  = NextOpenBrace( un_text, find );
    const size_t closebrace = MatchBrace( un_text, openbrace );

    const std::string suptext( un_text.begin() + openbrace+1,
                               un_text.begin() + closebrace );
    ans += 0.5 * EstimateLatexWidth( suptext );
    un_text.erase( un_text.begin() + find, un_text.begin() + closebrace + 1 );
  }

  // Dealing with super scripts and sub_scripts
  while( ( find = un_text.find( "_{" ) ) != std::string::npos ){
    const size_t openbrace  = NextOpenBrace( un_text, find );
    const size_t closebrace = MatchBrace( un_text, openbrace );

    const std::string subtext( un_text.begin() + openbrace+1,
                               un_text.begin() + closebrace );
    ans += 0.5 * EstimateLatexWidth( subtext );
    un_text.erase( un_text.begin() + find, un_text.begin() + closebrace + 1 );
  }

  // Removing left right decorators
  while( ( find = un_text.find( "#left" ) ) != std::string::npos ){
    un_text.erase( un_text.begin() + find, un_text.begin() + find + 5 );
  }

  while( ( find = un_text.find( "#right" ) ) != std::string::npos ){
    un_text.erase( un_text.begin() + find, un_text.begin() + find + 6 );
  }

  // Regex for finding generic decorators
  while( std::regex_match( un_text, match, dec_regex ) ){
    const std::string s = match[1];
    find = un_text.find( s );
    const size_t open  = NextOpenBrace( un_text, find );
    const size_t close = MatchBrace( un_text, open );
    un_text.erase( un_text.begin() + close );
    un_text.erase( un_text.begin() + find, un_text.begin() + s.length() );
  }

  // Replacing all symbols with just EM space
  while( std::regex_match( un_text, match, sym_regex ) ){
    std::regex_replace( un_text, sym_regex, "M" );
  }

  // Should have no more symbols
  for( const auto c : un_text ){
    if( char_width.find( c ) != char_width.end() ){
      ans += char_width.at( c );
    } else {
      ans += char_width.at( 'M' );// Adding M with for undetermined symbols
    }
  }

  return ans;
}

/**
 * @brief Estimate of the height of a Latex string in 'em' units.
 *
 * @details A line with no special characters is expected to be a line height of
 * ~1.08em (using the height ot the tallest regular character in ROOT: The
 * asterisk). A few processes are taken for fractions, for which we assume double
 * line height with some margin. One-line decorators (subscripts, supscripts and
 * square roots, simply adds an offset. The only special character that will the
 * taken care of would be the integral character (`\int`), which we assume to be
 * 2em tall.
 */
double
EstimateLatexHeight( const std::string& text )
{
  double ans = 1.08024691;// Maximum height of regular characters (*)
  std::string un_text( text );
  size_t find = std::string::npos;

  // Dealing with fractions in text
  while( ( find = un_text.find( "#frac{" ) ) != std::string::npos ){
    const size_t openbrace1  = NextOpenBrace( un_text, find );
    const size_t closebrace1 = MatchBrace( un_text, openbrace1 );
    const size_t openbrace2  = NextOpenBrace( un_text, closebrace1 );
    const size_t closebrace2 = MatchBrace( un_text, openbrace2 );

    const std::string numstr( un_text.begin()+openbrace1+1,
                              un_text.begin() + closebrace1 );
    const std::string denstr( un_text.begin()+openbrace2+1,
                              un_text.begin() + closebrace2 );
    ans = std::max( ans,
      EstimateLatexHeight( numstr )
      + EstimateLatexHeight( denstr )
      + fraclineheight );

    un_text.erase( un_text.begin() + find, un_text.begin() + closebrace2+1 );
  }

  // Dealing with integer symbols
  while( ( find = un_text.find( "#int" ) ) != std::string::npos ){
    ans = std::max( ans, int_height );
    un_text.erase( un_text.begin() + find, un_text.begin() + find + 4 );
  }

  // Only dealing with this for now....
  // Super scripts/ subscripts and parenthesis are nearly all covered in
  // the maximum height restriction. We are assuming no super complicated
  // Calculation to appear on plots.
  if( un_text.find( "#sqrt" ) != std::string::npos ){ ans += 0.22; }
  if( un_text.find( "^{" ) != std::string::npos ){ ans += 0.22; }
  if( un_text.find( "_{" ) != std::string::npos ){ ans += 0.22; }

  return ans;
}

}

}
