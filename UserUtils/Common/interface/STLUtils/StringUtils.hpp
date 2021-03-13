/**
 * @file    StringUtils.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Functions for string operations
 */
#ifndef USERUTILS_COMMON_STLUTILS_STRINGPATHUTILS_HPP
#define USERUTILS_COMMON_STLUTILS_STRINGPATHUTILS_HPP

#include <boost/format.hpp>
#include <string>
#include <vector>

namespace usr {

/**
 * @defgroup StringUtils String utilities
 * @brief    Helper string operation functions
 * @ingroup Common
 * @{
 */

/*-----------------------------------------------------------------------------
 * Getting the environment variable
   --------------------------------------------------------------------------*/
extern std::string GetEnv( const std::string& );

/*-----------------------------------------------------------------------------
 *  Globing related functions
   --------------------------------------------------------------------------*/
extern std::string GlobToRegex( const std::string& x );

/*-----------------------------------------------------------------------------
 *  Globbing in std::string format
   --------------------------------------------------------------------------*/
extern std::vector<std::string> GlobLocalStr( const std::string& pattern );

/*-----------------------------------------------------------------------------
 *  Generating random strings
   --------------------------------------------------------------------------*/
extern std::string RandomString( const unsigned length );

/*-----------------------------------------------------------------------------
 *  String operation functions
   --------------------------------------------------------------------------*/
extern void ToUpper( std::string& );
extern void StripSubstring( std::string&, const std::string& x );

extern size_t MatchBrace( const std::string&, const unsigned open_brace );
extern size_t NextOpenBrace( const std::string&, const unsigned start );


/**
 * @{
 * @brief Variadic interface for generating formmatted strings.
 *
 * The use of this function function would be more similar to the printf of C.
 * With the function returning a std::string instead of printing to a stream.
 * This is intended to simply boost::format calls were complicated formatting
 * isn't required. For examples:
 *
 * ```c++
 * std::cout << fstr("%d %d", myint, myunsigned ) << std::endl;
 * ```
 */
inline std::string
fstr( boost::format& fmt )
{  return boost::str( fmt ); }

template<typename TYPE,typename ... ARGS>
inline std::string
fstr( boost::format& fmt, TYPE&& T,  ARGS&& ... args )
{
  return fstr( fmt % std::forward<TYPE>(T), std::forward<ARGS>( args )... );
}

template<typename ... ARGS>
inline std::string
fstr( const std::string& fmt, ARGS&& ... args )
{
  boost::format f( fmt );
  return fstr( f, std::forward<ARGS>(args)... );
}
/** @} */

/** @} */

}/* usr */


#endif/* end of include guard: USERUTILS_COMMON_STLUTILS_STRINGPATHUTILS_HPP */
