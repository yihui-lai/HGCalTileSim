/**
 * @file RooArgContainer.hpp
 * @author Yi-Mu "Enoch" Chen
 * @brief RooArgContainter declaration class
 */

#ifndef USERUTILS_COMMON_ROOTUTILS_ROOARGCONTAINER_HPP
#define USERUTILS_COMMON_ROOTUTILS_ROOARGCONTAINER_HPP

#include <string>
#include <vector>

#include <RooCmdArg.h>

namespace usr
{

/**
 * @brief Simple container for helping with RooCmdArg parsing.
 */
class RooArgContainer : public std::vector<RooCmdArg>
{
public:
  static std::vector<std::string> CustomCommandList;
  static unsigned RegistorCommand( const std::string& );

  RooArgContainer( const std::vector<RooCmdArg>& arg_list,
                   const std::vector<RooCmdArg>& default_list = {} );
  virtual ~RooArgContainer();

  /**
   * @brief Template function for casting the result of Get() into a specific
   * type of RooCmdArg.
   *
   * This allow for access to more human readable version of the values stored in
   * the RooCmdArg.
   */
  template<typename TYPE>
  inline const TYPE
  Get() const { return TYPE( &Get( TYPE::CmdName ) ); }

  /**
   * @brief Template function for checking if a certain type of RooCmdArg exists
   * in the container. Reduces verbosity.
   */
  template<typename TYPE>
  inline bool
  Has() const { return Has( TYPE::CmdName ); }

  const RooCmdArg& Get( const std::string& name ) const;
  bool             Has( const std::string& name ) const;
  static bool      CheckList( const std::vector<RooCmdArg>&,
                             const std::string& name );
  RooLinkedList MakeRooList( const std::vector<std::string>& exclude = {} )
  const;
};

}// namespace

// MACRO FOR DEFINING A NEW COMMAND
#define USRUTILS_COMMON_REGISTERCMD( TYPE )              \
  const std::string TYPE::CmdName = #TYPE;               \
  static int __usrutils_common_register__dummy ## TYPE = \
     usr::RooArgContainer::RegistorCommand( #TYPE );

#endif
