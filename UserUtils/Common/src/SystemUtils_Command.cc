/**
 * @file    SystemUtils_Command.cc
 * @brief   Implementation of command line interaction functions
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/interface/SystemUtils/Command.hpp"
#include "UserUtils/Common/interface/SystemUtils/Time.hpp"
#else
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/SystemUtils/Command.hpp"
#include "UserUtils/Common/SystemUtils/Time.hpp"
#endif

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>

using namespace std;

namespace usr {

/**
 * @brief Returning the number of thread available on the machine. Simpler
 *        interface for the std::thread.
 */
unsigned
NumOfThreads()
{
  return std::thread::hardware_concurrency();
}

/**
 * @brief Executing a command, and storing the STD output to the return string.
 *
 * Notice that during the execution of this function, none of the outputs would
 * be displayed on screen, so the user is responsible for ensuring the command
 * would complete execution without user intervention.
 *
 * The main reference for the command could be obtained [here]( https://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c-using-posix)
 */
std::string
GetCMDOutput( const std::string& cmd )
{
  static const unsigned MAXBUFFER = 65536;
  static unsigned callcount       = 0;
  // Main reference
  std::array<char, MAXBUFFER> buffer;
  std::string result;
  std::shared_ptr<FILE> pipe( popen( cmd.c_str(), "r" ), pclose );
  ++callcount;
  if( !pipe ){ throw std::runtime_error( "popen() failed!" );}

  while( !feof( pipe.get() ) ){
    if( fgets( buffer.data(), 128, pipe.get() ) != nullptr ){
      result += buffer.data();
    }
  }

  return result;
}

/**
 * @brief The same as the GetCMDOutput() function, except that the std error
 *        outputs would be discarded into `/dev/null`.
 */
std::string
GetCMDSTDOutput( const std::string& cmd )
{
  return GetCMDOutput( cmd+" 2> /dev/null" );

}

/**
 * @brief The same as the GetCMDOutput() function, except that only the
 *        std error output would be stored, and the std output would be
 *        discarded  into `/dev/null`.
 */
std::string
GetCMDERROutput( const string& cmd )
{
  return GetCMDOutput( cmd+" > /dev/null  2>&1 " );
}

/**
 * @brief Return the number of the processes whose name contains the string
 *        "cmd", but doesn't include the string "exclude".
 *
 * This function uses standard UNIX tools `ps`, `awk`, `wc` and `grep` to
 * get the desired results.
 */
int
HasProcess( const std::string& cmd, const std::string& exclude )
{
  static const std::string common = "ps -U ${USER} | awk '{print $4}' ";
  static const std::string count  = " | wc --lines";
  const std::string grepcmd       = " | grep " + cmd;
  const std::string excludecmd    = exclude == "" ? "" :
                                    " | grep --invert-match " + exclude;

  const std::string ans = GetCMDOutput( common + grepcmd + excludecmd + count );
  return stoi( ans );
}

/**
 * @brief Sleep current thread until the results of HasProcess returns 0
 *
 * To avoid confusion during the waiting, this functions is made intensionally
 * noisy, printing the machine, system time, and number of instances of the
 * specific command each time the HasProcess is checked for a zero return.
 * By default, the check is performed every 100ms.
 */
void
WaitProcess(
  const std::string& x,
  const std::string& exclude,
  const unsigned     sleeptime )
{
  while( 1 ){
    const int proc_count          = HasProcess( x, exclude );
    const std::string time_string = CurrentTime();

    if( proc_count == 0 ){ break; }

    cout <<
      usr::fstr( "\r[%s|%s] %d instance(s) of [%s]...",
      GetEnv( "HOSTNAME" ),
      time_string,
      proc_count,
      x )
         << flush;

    SleepMillSec( sleeptime );
  }

  cout << "All Done!" << endl;
}

}/* usr */
