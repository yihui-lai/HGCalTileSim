/**
 * @file    ArgumentExtender.cc
 * @brief   Implementations for the ArgumentExtender class
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/ArgumentExtender.hpp"
#include "UserUtils/Common/interface/Format.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#else
#include "UserUtils/Common/ArgumentExtender.hpp"
#include "UserUtils/Common/Format.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#endif

#include <boost/algorithm/string.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options/errors.hpp>

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iostream>

namespace pt  = boost::property_tree;
namespace opt = boost::program_options;
namespace fs  = std::experimental::filesystem;

namespace usr {

/**
 * @brief Construction of the class property tree with a list of files.
 *
 * The constructor will ensure that the json structure for each entry under a
 * given options is uniformed in the json file. Notice that the options listed
 * in the json file will be considered a required arguments (hence the support
 * for multiple json files).
 *
 * An exception would be thrown if any parsing error occurs. Multiple options
 * with the same name would also be considered an exception.
 *
 * After successfully parsing the json file. The option `help/h` would be added
 * to for use (if used, the class will diplay the defined options and exit).
 */
void
ArgumentExtender::_init( const std::vector<std::string>& filelist )
{
  _exttree = FromJsonFiles( filelist );

  // Parsing the format of each options listed in ext trees
  BOOST_FOREACH( const pt::ptree::value_type& child, Tree() ){
    const std::string optname = child.first.data();
    if( child.second.size() == 0 ){
      std::cerr << "Bad format for option ["<< optname << "]" << std::endl;
      throw std::invalid_argument( optname );
    }

    std::vector<std::string> exttaglist;

    BOOST_FOREACH( const pt::ptree::value_type& subchild, child.second ){
      const std::string optval = subchild.first.data();
      if( subchild.second.size() == 0 ){
        std::cerr << "Bad format for option ["<< optname << "] "
                  << "with value [" << optval << "]" << std::endl;
        throw std::invalid_argument( optval );
      }

      std::vector<std::string> this_exttaglist;
      BOOST_FOREACH( const pt::ptree::value_type& leafp, subchild.second ){
        this_exttaglist.push_back( leafp.first.data() );
      }
      std::sort( this_exttaglist.begin(), this_exttaglist.end() );
      if( this_exttaglist.size() == 0 ){
        std::cerr << "No extended values for option ["<< optname << "] "
                  << "with value [" << optval << "]" << std::endl;
        throw std::invalid_argument( optval );
      }

      if( exttaglist.size() == 0  ){// Do not compare if this is empty
        exttaglist = this_exttaglist;
      } else {

        if( exttaglist.size() > this_exttaglist.size() ){
          for( const auto& tag : exttaglist ){
            if( !FindValue( this_exttaglist, tag ) ){
              std::cerr << "Missing extended value [" << tag <<"] "
                        << "for option ["<< optname << "] "
                        << "with value [" << optval << "]" << std::endl;
              throw std::invalid_argument( tag );
            }
          }
        } else if( exttaglist.size() < this_exttaglist.size() ){
          for( const auto& tag : this_exttaglist ){
            if( !FindValue( this_exttaglist, tag ) ){
              std::cerr << "Missing extended value [" << tag <<"] "
                        << "for option ["<< optname << "] "
                        << "with value [" << optval << "]" << std::endl;
              throw std::invalid_argument( tag );
            }
          }
        } else {
          for( size_t i = 0; i < exttaglist.size(); ++i ){
            if( exttaglist.at( i ) != this_exttaglist.at( i ) ){
              const std::string misval = exttaglist.at( i );
              std::cerr << "Undefined extended value [" << misval <<"] "
                        << "for option ["<< optname << "] "
                        << "with value [" << optval << "]" << std::endl;
              throw std::invalid_argument( misval );
            }
          }
        }
      }
    }
  }

  // Adding default arguments.
  Description().add_options()
    ( "help,h", "print help options and exit program" )
  ;
}

/**
 * @brief default destructor (nothing, relying on member auto cleaning)
 */
ArgumentExtender::~ArgumentExtender(){}

/**
 * @brief Adding more options definition to the class.
 *
 * Currently the only way of adding more options is to manually define
 * a new options_description instance, and adding that to the class.
 * The internal options_description instance would only directly contain
 * the options definition of "help/h".
 */
ArgumentExtender&
ArgumentExtender::AddOptions( const opt::options_description& desc )
{
  Description().add( desc );
  return *this;
}

// Parse the program input given the standard argc, argv interface.
// Throw exception when parse error or exit the program if "help" options
// is received.
// All functions listed under the json files would be listed as mandetory.
// and a parse error would be thrown if input is not detected.

/**
 * @brief parsing the program input with the standard argc, argv interface.
 *
 * Raising exception if there is any error in parsing. Additional criteria
 * that counts as an error in parsing includes:
 * 1. The json file contains an options not listed in the final
 *     options_description instance.
 * 2. A option listed in the json does not have a corresponding user input.
 * 3. The user input for an option listed in the json file doesn't match
 *    any of inputs listed in the json file.
 *
 * After all parsing is complete, the function will check for the "help/h"
 * input. If it exists, the function will print the options_description and
 * terminate the program.
 */
void
ArgumentExtender::ParseOptions( int argc, char** argv )
{
  try {
    opt::store( opt::parse_command_line( argc, argv, Description() ), Args() );
    opt::notify( Args() );
  } catch( boost::exception& e ){
    std::cerr << "Error parsing command!" << std::endl
              << boost::diagnostic_information( e )
              << Description()
              << std::endl;
    throw;// Continue throwing exception
  }

  // Checking that each of the options listed in the json files
  // are present. Will need to
  BOOST_FOREACH( const pt::ptree::value_type& child, Tree() ){
    const std::string optname = child.first.data();

    try {
      Description().find( optname, true );
    } catch( boost::exception& e ){
      std::cerr << "Options [" << optname << "] defined in file"
                << "does not have description defined!" << std::endl;
      throw;
    }

    if( !CheckArg( optname ) ){
      std::cerr << "Missing options for [" << optname << "]" << std::endl
                << Description()
                << std::endl;
      BOOST_THROW_EXCEPTION( po::required_option( optname ) );
    }

    std::vector<std::string> validoptionlist;
    if( !CheckQuery( Tree(), optname, Arg<std::string>( optname ) ) ){
      std::cerr << "Extended values for options for [" << optname << "] "
                << " with value [" << Arg<std::string>( optname ) << "] "
                << " is not defined!"
                << std::endl
                << "Available values: " << std::flush;
      BOOST_FOREACH( const pt::ptree::value_type& subchild, child.second ){
        const std::string optvalue = subchild.first.data();
        std::cerr << "[" << optvalue << "] " << std::flush;
      }
      std::cerr << std::endl;

      throw std::invalid_argument( Arg<std::string>( optname ) );
    }
  }

  if( CheckArg( "help" ) ){
    PrintHelpAndExit();
  }
}

/**
 * @brief Printing help message of arguments and exit the program.
 */
void
ArgumentExtender::PrintHelpAndExit() const
{
  std::cerr << Description() << std::endl;
  std::exit( EXIT_SUCCESS );
}

/**
 * @brief Checking if the user input exists for an option.
 * @details Should be called after the ParseOptions() method has been called.
 */
bool
ArgumentExtender::CheckArg( const std::string& opt ) const
{
  return Args().count( opt );
}

/**
 * @brief Setting the file prefix for the file path generation
 *
 * This is the directory path that will be used regardless of any input status.
 */
void
ArgumentExtender::SetFilePrefix( const fs::path pre )
{
  _prefix = pre;
}

/**
 * @brief Setting the directory naming scheme list
 *
 * This will **remove** any existing scheme for the directory naming.
 */
void
ArgumentExtender::SetDirScheme( const PathScheme& newscheme )
{
  _dirscheme = newscheme;
}

/**
 * @brief Incrementally adding single directory name scheme to the final list
 */
void
ArgumentExtender::AddDirScheme( const ArgPathScheme& arg )
{
  _dirscheme.push_back( arg );
}

/**
 * @brief Adding a list of directory naming scheme to the existing list.
 *
 * This is extend the existing list instead of straight up replacing it.
 */
void
ArgumentExtender::AddDirScheme( const PathScheme& newscheme )
{
  for( const auto& x : newscheme ){
    _dirscheme.push_back( x );
  }
}

/**
 * @brief Setting the basename naming scheme list
 *
 * This will **remove** any existing scheme for the directory naming.
 */
void
ArgumentExtender::SetNameScheme( const PathScheme& newscheme )
{
  _namescheme = newscheme;
}

/**
 * @brief Incrementally adding single basename name scheme to the final list
 */
void
ArgumentExtender::AddNameScheme( const ArgPathScheme& arg )
{
  _namescheme.push_back( arg );
}

/**
 * @brief Adding a list of basename naming scheme to the existing list.
 *
 * This is extend the existing list instead of straight up replacing it.
 */
void
ArgumentExtender::AddNameScheme( const PathScheme& newscheme )
{
  for( const auto& x : newscheme ){
    _namescheme.push_back( x );
  }
}

/**
 * @brief Generating a path based of existing naming scheme in the class
 *
 * The directory scheme is to use different inputs to generated nested directory
 * names, while the name scheme will generate the resulting base name of the
 * file. The character separating the option inputs values will simply be '_'.
 *
 * For details on the generation of individual options string, see the private
 * method `genPathString()`.
 */
fs::path
ArgumentExtender::MakeFile(
  const std::string& nameprefix,
  const std::string& ext ) const
{
  fs::path dirname  = _prefix;
  fs::path filename = nameprefix;

  for( const auto& x : _dirscheme ){
    dirname /= genPathString( x );
  }

  for( const auto& x : _namescheme ){
    const std::string str = genPathString( x );
    if( str != "" ){
      filename += "_" + genPathString( x );
    }
  }

  filename += "." + ext;

  return fs::path( dirname / filename );
}

fs::path
ArgumentExtender::MakePDFFile( const std::string& x ) const
{
  return MakeFile( x, "pdf" );
}

fs::path
ArgumentExtender::MakePNGFile( const std::string& x ) const
{
  return MakeFile( x, "png" );
}

fs::path
ArgumentExtender::MakeTXTFile( const std::string& x ) const
{
  return MakeFile( x, "txt" );
}

fs::path
ArgumentExtender::MakeTEXFile( const std::string& x ) const
{
  return MakeFile( x, "tex" );
}

/**
 * @brief Generating the string for a single option input.
 *
 * For the given scheme of { option, prefixstring }, the function will return the
 * string: prefixstring + option_input_as_string. Some substitutions will be name
 * to make the output more friendly to command lines:
 *
 * 1. Decimal points will be replaced with 'p'
 * 2. Spaces (for list multitoken arguments) will be replaced with '-'
 * 3. If the option doesn't have a input (a.k.a. The option is just a boolean
 *    flag), the the option_input_as_string is would be "On/Off" depending on the
 *    flag. If prefixstring if empty, then the return string would be "option" or
 *    ""  depending on the options.
 */
std::string
ArgumentExtender::genPathString( const ArgPathScheme& x ) const
{
  std::string inputstring = "";

  if( IsBooleanFlag( x.option ) ){
    inputstring = genPathString_Boolean( x );
  } else if( IsMultiToken( x.option ) ){
    inputstring = genPathString_List( x );
  } else {
    inputstring = genPathString_Single( x );
  }

  boost::replace_all( inputstring, ".", "p" );
  boost::replace_all( inputstring, " ", "-" );

  return inputstring;
}

std::string
ArgumentExtender::genPathString_Boolean( const ArgPathScheme& x ) const
{
  const bool flag = CheckArg( x.option );

  if( x.pathstring == "" ){
    return flag ? x.option : "";
  } else {
    return x.pathstring + std::string( flag ? "On" : "Off" );
  }
}

std::string
ArgumentExtender::genPathString_List( const ArgPathScheme& x ) const
{
  if( !CheckArg( x.option ) ){ return ""; }

  std::string ans = "";

  if( Args()[x.option].value().type() == typeid( std::vector<std::string> ) ){
    for( const auto in : ArgList<std::string>( x.option ) ){
      ans += in + '-';
    }

    ans.erase( ans.end() - 1 );
  } else if( Args()[x.option].value().type() == typeid( std::vector<int> ) ){
    for( const auto in : ArgList<int>( x.option ) ){
      ans += std::to_string( in ) + '-';
    }

    ans.erase( ans.end() - 1 );
  } else if( Args()[x.option].value().type() == typeid( std::vector<double> ) ){
    for( const auto in : ArgList<double>( x.option ) ){
      ans += usr::fmt::base::decimal( in ).str() + '-';
    }

    ans.erase( ans.end() - 1 );
  }

  return x.pathstring + ans;
}

std::string
ArgumentExtender::genPathString_Single( const ArgPathScheme& x ) const
{
  if( !CheckArg( x.option ) ){ return ""; }

  std::string ans = "";

  if( Args()[x.option].value().type() == typeid( std::string ) ){
    ans = Arg<std::string>( x.option );
  } else if( Args()[x.option].value().type() == typeid( int ) ){
    ans = std::to_string( Arg<int>(x.option) );
  } else if( Args()[x.option].value().type() == typeid( double ) ){
    ans = usr::fmt::base::decimal( Arg<double>(x.option) ).str();
  }

  return x.pathstring + ans;
}


/**
 * @brief Checking if a certain options is a bool flag option (i.e. declared with
 *        no semantic information)
 *
 * Options declare with no semantic information (i.e. no
 * boost::program_options::value<std::string>()] and the like), are typically
 * used as boolean flags. The function will return if a certain options is
 * declared as such, using the semantic->mix_token and max_token method (both are
 * zero if it is)
 *
 * In the case that the option doesn't exist in the class instance, false is
 * always returned.
 */
bool
ArgumentExtender::IsBooleanFlag( const std::string& x ) const
{
  const auto description = Description().find_nothrow( x, true );

  if( !description ){ return false; }// Not found

  return description->semantic()->min_tokens() == 0 &&
         description->semantic()->max_tokens() == 0;
}

/**
 * @brief Checking if a certain options takes in multiple inputs.
 *
 * Options declare with the multi-token semantic information with have a non-1
 * max_token method in the sematic information
 *
 * In the case that the option doesn't exist in the class instance, false is
 * always returned.
 */
bool
ArgumentExtender::IsMultiToken( const std::string& x ) const
{
  const auto description = Description().find_nothrow( x, true );

  if( !description ){ return false; }// Not found

  return description->semantic()->max_tokens() > 1;
}


}/* usr */
