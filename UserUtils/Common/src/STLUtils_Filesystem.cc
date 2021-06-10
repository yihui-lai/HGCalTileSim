/**
 * @file    STLUtils_Filesystem.cc
 * @brief   Implementation of filesystem enhancement functions
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/Filesystem.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#else
#include "UserUtils/Common/STLUtils/Filesystem.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#endif

#include <glob.h>
#include <iostream>

namespace usr  {

/**
 * @brief returning path to the CMSSW/src directory.
 */
fs::path
cmssw_src(){ return GetEnv( "CMSSW_BASE" ) + "/src/";  }

/**
 * @brief returning path to a main package directory.
 */
fs::path
pkgpath( const std::string& x )
{
  return cmssw_src() / fs::path( x );
}

/**
 * @brief returning path to a sub-package directory
 */
fs::path
subpkgpath( const std::string& x, const std::string& y )
{
  return pkgpath( x )/fs::path( y );
}

/**
 * @brief returning path to the data directory of a sub-package.
 */
fs::path
datapath( const std::string& x, const std::string& y )
{
  return subpkgpath( x, y )/fs::path( "data" );
}

/**
 * @brief returning path to the results directory of a sub-package (unofficial)
 */
fs::path
resultpath( const std::string& x, const std::string& y )
{
  return subpkgpath( x, y )/fs::path( "result" );
}

/**
 * @brief returning the globbing results of a query string.
 *
 * Currently using the POSIX glob.h library, which is always available for the
 * CMSSW environment. This functions is simple used for translating to the
 * std::filesystem format.
 */
std::vector<fs::path>
GlobLocal( const std::string& query )
{
  glob_t glob_result;
  glob( query.c_str(), GLOB_TILDE, NULL, &glob_result );

  std::vector<fs::path> ans(
    glob_result.gl_pathv,
    glob_result.gl_pathv + glob_result.gl_pathc
    );
  globfree( &glob_result );

  return ans;
}


/**
 * @brief Checking if the parent directory of the filepath exists, and
 *        attempting to create if it doesn't exists.
 *
 * This functions is intentionally noisy when creating the parent directories
 * as this calls system commands and is potentially dangerous.
 */
bool
MakeParent( const fs::path& filepath )
{
  const fs::path parent = fs::path( filepath ).parent_path();
  if( fs::is_directory( parent ) || parent.empty() ){
    return true;// Early exit if already exists
  }

  if( fs::create_directories( parent ) ){
    std::cout << "New directory [" << parent << "] created!"
              << std::endl;
    return true;
  } else {
    std::cout << "Failed making new directory [" << parent << "]!"
              <<  std::endl;
    return false;
  }
}

}/* usr  */
