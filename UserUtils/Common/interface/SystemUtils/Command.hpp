/**
 * @file    Command.hpp
 * @brief   Functions for command-line interaction
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#ifndef USERUTILS_COMMON_SYSTEMUTILS_COMMAND_HPP
#define USERUTILS_COMMON_SYSTEMUTILS_COMMAND_HPP

#include <string>

namespace usr {

/**
 * @defgroup SystemCommands System Commands
 * @brief    Helper functions for interacting with system commands.
 * @ingroup  Common
 * @{
 */

/*-----------------------------------------------------------------------------
 *  General helper functions
   --------------------------------------------------------------------------*/
extern unsigned NumOfThreads();

/*-----------------------------------------------------------------------------
 *  Execute command and output results.
   --------------------------------------------------------------------------*/
extern std::string GetCMDOutput( const std::string& );
extern std::string GetCMDSTDOutput( const std::string& );
extern std::string GetCMDERROutput( const std::string& );

/*-----------------------------------------------------------------------------
 *  Existing process detection and waiting.
   --------------------------------------------------------------------------*/
extern int HasProcess( const std::string&, const std::string& = "" );

extern void WaitProcess(
  const std::string& process_tag,
  const std::string& exclude = "",
  const unsigned     waitmills = 100
  );
/* @} */


} /* usr */


#endif /* end of include guard: USERUTILS_COMMON_SYSTEMUTILS_COMMAND_HPP */
