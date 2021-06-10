/**
 * @file   BoostUtils_PTreeUtils.cc
 * @author [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief  Implementing boost property tree operations
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/BoostUtils/PTreeUtils.hpp"
#else
#include "UserUtils/Common/BoostUtils/PTreeUtils.hpp"
#endif

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <queue>

namespace usr {

/**
 * @brief populating a property tree from a json file.
 */
pt::ptree
FromJsonFile( const std::string& filename )
{
  pt::ptree ans;
  pt::read_json( filename, ans );
  return ans;
}

/**
 * @brief populating a property tree from a list of json files.
 * @details Calls the MergeTrees() functions for tree created by the different
 * files.
 */
pt::ptree
FromJsonFiles( const std::vector<std::string>& filelist )
{
  pt::ptree ans;

  for( const auto& file : filelist ){
    pt::ptree tmp;
    pt::read_json( file, tmp );
    ans = MergeTrees( ans, tmp );
  }

  return ans;
}

/**
 * @brief merging the contents of two property trees
 *
 * Main reference from [here](https://paste.tbee-clan.de/TX2Vm).
 * The second tree is traversed using BFS algorithm and have it's contents
 * pushed accordingly. Conflicting entries may have undesirable behaviour.
 * Currently untested.
 */
pt::ptree
MergeTrees( const pt::ptree& first, const pt::ptree& second )
{
  pt::ptree merged = first;

  std::queue<pt::ptree> qvalues;
  qvalues.push( second );
  std::queue<std::string> qkeys;

  while( !qvalues.empty() ){
    pt::ptree current = qvalues.front();
    qvalues.pop();
    std::string keychain = "";
    if( !qkeys.empty() ){
      keychain = qkeys.front();
      qkeys.pop();
    }

    BOOST_FOREACH( const pt::ptree::value_type& child, current ){
      if( child.second.size() == 0 ){// Leaf
        // No "." for first level entries
        const std::string s =
          keychain == "" ? child.first.data() :
          keychain + "." + child.first.data();
        merged.put( s, child.second.data() );
      } else {// Subtree
        // Put keys (identifiers of subtrees) and all of its parents
        // aside for later iteration. Keys on first level have no parents
        const std::string newkey
          = keychain == "" ? child.first.data() :
            keychain + "." + child.first.data();
        qkeys.push( newkey );
        qvalues.push( child.second );
      }
    }// -- End of BOOST_FOREACH
  }// --- End of while

  return merged;
}

/**
 * @brief printing a property tree using recursive algorithms. The level
 *        automatically increments by 1 each level into the tree.
 *
 * The user can specify a level greater than 0, but the only result is the
 * the screen output is shifted to the right.
 */
void
PrintPTree( const pt::ptree& tree, unsigned level )
{
  for( const auto& it : tree ){
    for( unsigned i = 0; i < level; ++i ){
      std::cout << "\t" << std::flush;
    }

    std::cout << it.first << ": "
              << it.second.get_value<std::string>() << std::endl;
    PrintPTree( it.second, level + 1 );
  }
}

/**
 * @brief creating a query string based on a list of sub-level tags.
 *
 * Given a list of strings {string1, string2, string3}, the return
 * string will concat the strings into a single property tree query as:
 * "string1.string2.string3".
 */
std::string
MakeQueryString( const std::vector<std::string>& list )
{
  return boost::join( list, "." );
}

/**
 * @brief checking if anything exists in the property tree at the query point.
 */
bool
CheckQuery( const pt::ptree& tree, const std::string& query )
{
  try {
    tree.get_child(query);
    return true;
  } catch( boost::exception& e ){
    return false;
  }
}

/**
 * @brief returning the sub-tree at the query point.
 */
pt::ptree
GetSubTree( const pt::ptree& tree, const std::string& query )
{
  return tree.get_child( query );
}

}/* usr */
