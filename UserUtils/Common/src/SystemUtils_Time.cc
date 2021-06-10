/**
 * @file    SystemUtils_Time.cc
 * @brief   Implementation of system time related functions
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#include <chrono>
#include <ctime>
#include <string>
#include <thread>

namespace usr {

/**
 * @brief unix time template functions for helping with the various time units
 *        types.
 */
template<typename T>
inline size_t
CurrentTimeIn()
{
  return std::chrono::duration_cast<T>(
    std::chrono::high_resolution_clock::now().time_since_epoch()
    ).count();
}

/**
 * @brief Current unix time in nano seconds.
 */
size_t
CurrentTimeInNanSec()
{
  return CurrentTimeIn<std::chrono::nanoseconds>();
}

/**
 * @brief current unix time in microseconds
 */
size_t
CurrentTimeInMuSec()
{
  return CurrentTimeIn<std::chrono::microseconds>();
}

/**
 * @brief Current unix time in milliseconds
 */
size_t
CurrentTimeInMSec()
{
  return CurrentTimeIn<std::chrono::milliseconds>();
}

/**
 * @brief Current unix time in seconds.
 */
size_t
CurrentTimeInSec()
{
  return CurrentTimeIn<std::chrono::seconds>();
}

/**
 * @brief Representing the current time as as string.
 *
 * Currently there is still not official C++ support, so this function uses
 * the old C API. Details on how the `fmt` string represents the format could
 * be found [here](http://www.cplusplus.com/reference/ctime/strftime/)
 */
std::string
CurrentTime( const std::string& fmt )
{
  static const unsigned MAXBUFFER = 80;
  time_t rawtime;
  struct tm* timeinfo;
  char buffer [MAXBUFFER];
  time( &rawtime );
  timeinfo = localtime( &rawtime );
  strftime( buffer, MAXBUFFER, fmt.c_str(), timeinfo );
  return buffer;
}


/**
 * @brief Sleeping template function for a handling various time types.
 */
template<typename T>
inline void
SleepFor( const unsigned x )
{
  std::this_thread::sleep_for( T( x ) );
}

/**
 * @brief sleeping thread for x microseconds.
 */
void
SleepMicroSec( const unsigned x )
{
  SleepFor<std::chrono::microseconds>( x );
}

/**
 * @brief sleeping thread for x milliseconds.
 */
void
SleepMillSec( const unsigned x )
{
  SleepFor<std::chrono::milliseconds>( x );
}

/**
 * @brief sleeping thread for x seconds.
 */
void
SleepSec( const unsigned x )
{
  SleepFor<std::chrono::seconds>( x );
}

}/* usr */
