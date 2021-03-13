/**
 * @file    STLUtils_OStreamUtils.cc
 * @brief   ostream interaction class implementation.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/OStreamUtils.hpp"
#else
#include "UserUtils/Common/STLUtils/OStreamUtils.hpp"
#endif

#include <ostream>

namespace usr {

/**
 * @details currently the clear line implementation consists of printing
 *          multiple backspace and whitespaces. Maybe not the best?
 */
std::ostream&
operator<<( std::ostream& os, const clearline& )
{
  for( size_t i = 0; i < 255; ++i ){
    os << "\b \b" << std::flush;
  }

  return os;
}

separator::separator( const char x, const size_t m ) :
  token( x ),
  n(m){}

separator::~separator(){}

/**
 * @brief Implementation is essentially what it says on the tin
 *
 * duplicating the character multiple times.
 */
std::ostream&
operator<<( std::ostream& os, const separator& x )
{
  for( size_t i = 0; i < x.n; ++i ){
    os << x.token << std::flush;
  }

  return os;
}

}/* usr */
