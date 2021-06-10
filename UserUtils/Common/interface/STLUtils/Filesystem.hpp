/**
 * @file    Filesystem.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Defining functions for helping with file system interaction
 */
#ifndef USERUTILS_COMMON_STLUTILS_FILESYSTEM_HPP
#define USERUTILS_COMMON_STLUTILS_FILESYSTEM_HPP

#include <string>
#include <vector>
#include <experimental/filesystem>

namespace usr {

namespace fs=std::experimental::filesystem;

/**
 * @defgroup filesys Filesystem helper
 * @brief Simplification for the std::filesystem data structure.
 * @ingroup Common
 * @{
 */
/*-----------------------------------------------------------------------------
 *  Commonly used paths.
   --------------------------------------------------------------------------*/
extern fs::path cmssw_src();
extern fs::path pkgpath( const std::string& );
extern fs::path subpkgpath( const std::string&, const std::string& );
extern fs::path datapath( const std::string&, const std::string& );
extern fs::path resultpath( const std::string&, const std::string& );

/*-----------------------------------------------------------------------------
 *  Globbing
   --------------------------------------------------------------------------*/
extern std::vector<fs::path> GlobLocal( const std::string& );

/*-----------------------------------------------------------------------------
 *  File checking functions
   --------------------------------------------------------------------------*/
extern bool MakeParent( const fs::path& );

/** @*/
}

#endif /* end of include guard: USERUTILS_COMMON_STLUTILS_FILESYSTEM_HPP */
