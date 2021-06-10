/**
 * @file    Time.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifndef USERUTILS_COMMON_SYSTEMUTILS_TIMP_HPP
#define USERUTILS_COMMON_SYSTEMUTILS_TIMP_HPP

#include <string>

namespace usr {

/**
 * @defgroup SystemTime System Time
 * @brief    Simple interface to the std::chrono for time handling.
 * @ingroup  Common
 * @{
 */
/*-----------------------------------------------------------------------------
 *  Getting UNIX time as a 64-bit number
   --------------------------------------------------------------------------*/
extern size_t CurrentTimeInNanSec();
extern size_t CurrentTimeInMuSec();
extern size_t CurrentTimeInMSec();
extern size_t CurrentTimeInSec();

/*-----------------------------------------------------------------------------
 *  Displaying present time
   --------------------------------------------------------------------------*/
extern std::string CurrentTime( const std::string& fmt="%Y-%m-%d %X");

/*-----------------------------------------------------------------------------
 *  Sleeping a thread for various seconds
   --------------------------------------------------------------------------*/
extern void SleepMicroSec( const unsigned );
extern void SleepMillSec( const unsigned );
extern void SleepSec( const unsigned );

/* @} */

} /* usr */

#endif/* end of include guard: USERUTILS_COMMON_SYSTEMUTILS_TIMP_HPP */
