/**
 * @file    PTreeUtils.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Defining common boost property tree operations
 */
#ifndef USERUTILS_COMMON_BOOSTUTILS_PTREEUTILS_HPP
#define USERUTILS_COMMON_BOOSTUTILS_PTREEUTILS_HPP

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#else
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#endif

#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace usr {

/**
 * @defgroup ptreeutils Boost properties tree helper
 * @brief Simple functions aimed to simplifying the use of the boost property
 *        tree
 * @ingroup Common
 * @{
 */
namespace pt = boost::property_tree;

/*-----------------------------------------------------------------------------
 *  Making ptree objects from a given file
   --------------------------------------------------------------------------*/
extern pt::ptree FromJsonFile( const std::string& );
extern pt::ptree FromJsonFiles( const std::vector<std::string>& );
// variadic interface
template<typename ... TS>
extern pt::ptree FromJsonFiles( const std::string& x, TS ... args );

/*-----------------------------------------------------------------------------
 *  Two way merging functions.
   --------------------------------------------------------------------------*/
extern pt::ptree MergeTrees( const pt::ptree& first, const pt::ptree& second );

/*-----------------------------------------------------------------------------
 *  Dumping the contents of a pt::ptree object
   --------------------------------------------------------------------------*/
extern void PrintPTree( const pt::ptree&, unsigned level = 0 );

/*-----------------------------------------------------------------------------
 *  Query functions
   --------------------------------------------------------------------------*/
extern std::string MakeQueryString( const std::vector<std::string>& );

extern bool CheckQuery( const pt::ptree&, const std::string& );

template<typename T>
extern T GetSingle( const pt::ptree& tree, const std::string& query );

template<typename T>
extern std::vector<T> GetList( const pt::ptree&, const std::string& );

extern pt::ptree GetSubTree( const pt::ptree&, const std::string& );

/*-----------------------------------------------------------------------------
 *  Variadic interface
   --------------------------------------------------------------------------*/
template<typename ... TS>
std::string MakeQueryString( const std::string&, TS ... args );

template<typename ... TS>
extern bool CheckQuery( const pt::ptree&, const std::string&, TS ... args );

template<typename T, typename ... TS>
extern T GetSingle( const pt::ptree&, const std::string&, TS ... args );

template<typename T, typename ... TS>
extern std::vector<T> GetList( const pt::ptree&, const std::string&, TS ... args );

template<typename ... TS>
extern pt::ptree GetSubTree( const pt::ptree&, const std::string&, TS ... args );

/** @} */

/*-----------------------------------------------------------------------------
 *  template implementations
   --------------------------------------------------------------------------*/

/**
 * @brief getting single data at given address.
 *
 * This assumes that the address points to a single data type and not a
 * sub-tree. Boost handles the conversion and exception raising.
 */
template<typename T>
extern T
GetSingle( const pt::ptree& tree, const std::string& query )
{
  return tree.get<T>( query );
}

/**
 * @brief Getting list of objects at given address.
 *
 * Main reference of code came from [here](http://stackoverflow.com/questions/23481262/using-boost-property-tree-to-read-int-array)
 */
template<typename T>
extern std::vector<T>
GetList( const pt::ptree& tree, const std::string& query )
{
  std::vector<T> ans;

  for( auto& v : tree.get_child( query ) ){
    if( !v.first.empty() ){
      std::cerr << "Warning! Skipping over illegal format at branch: (" << query
                << ")  with index value: (" << v.first.data() << ")" << std::endl;
      continue;
    }
    ans.push_back( v.second.get_value<T>() );
  }

  return ans;
}

/**
 * @brief variadic interface to allow construction of property tree from
 * arbitrary number of input files.
 */
template<typename ... TS>
extern pt::ptree
FromJsonFiles( const std::string& x, TS ... args )
{
  return FromJsonFiles( MakeVector<std::string>( x, args ... ) );
}

/**
 * @brief variadic function for the MakeQueryString() function.
 */
template<typename ... TS>
std::string
MakeQueryString(
  const std::string& x, TS ... args )
{
  return usr::MakeQueryString( usr::MakeVector<std::string>( x, args ... ) );
}


/**
 * @brief variadic function to the CheckQuery() function. all string arguments
 * would be used to create the single query string.
 */
template<typename ... TS>
extern bool
CheckQuery( const pt::ptree& tree, const std::string& x, TS ... args )
{
  return usr::CheckQuery( tree, usr::MakeQueryString( x, args ... ) );
}

/**
 * @brief variadic interface to the GetSingle() function. All string arguments
 * would be used to create the single query string.
 */
template<typename T, typename ... TS>
extern T
GetSingle( const pt::ptree& tree, const std::string& x, TS ... args )
{
  return usr::GetSingle<T>( tree, usr::MakeQueryString( x, args ... ) );
}

/**
 * @brief variadic interface to the GetList() function. All string arguments
 * would be used to create the single query string.
 */
template<typename T, typename ... TS>
extern std::vector<T>
GetList( const pt::ptree& tree, const std::string& x, TS ... args )
{
  return usr::GetList<T>( tree, usr::MakeQueryString( x, args ... ) );
}

/**
 * @brief variadic interface to the GetSubTree() function. All string arguments
 * would be used to create the single query string.
 */
template<typename ... TS>
extern pt::ptree
GetSubTree( const pt::ptree& tree, const std::string& x, TS ... args )
{
  return usr::GetSubTree( tree, usr::MakeQueryString( x, args ... ) );
}

}/* usr*/

#endif/* end of include guard: USERUTILS_COMMON_BOOSTUTILS_PTREEUTILS_HPP */
