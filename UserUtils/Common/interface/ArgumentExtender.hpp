/**
 * @file    ArgumentExtender.hpp
 * @brief   Class for managing main function inputs.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */
#ifndef USERUTILS_COMMON_ARGUMENTEXTENDER_HPP
#define USERUTILS_COMMON_ARGUMENTEXTENDER_HPP

#include <boost/program_options.hpp>
#include <experimental/filesystem>
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "UserUtils/Common/interface/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/interface/STLUtils/VectorUtils.hpp"
#else
#include "UserUtils/Common/BoostUtils/PTreeUtils.hpp"
#include "UserUtils/Common/STLUtils/StringUtils.hpp"
#include "UserUtils/Common/STLUtils/VectorUtils.hpp"
#endif

namespace usr {

namespace po = boost::program_options;
namespace fs = std::experimental::filesystem;

class ArgumentExtender
{
public:

  /**
   *@ brief Initialization without additional json file
   */
  ArgumentExtender(){ _init( {} ); }

  /**
   * @brief initialization from single json file.
   * @details see ArgumentExtender::_init() for more details.
   */
  ArgumentExtender( const std::string& filename )
  { _init( {filename} ); }

  /**
   * @brief initialization from multiple json files
   * @details see, ArgumentExtender::_init() for more details.
   */
  ArgumentExtender( const std::vector<std::string>& list )
  { _init( list ); }

  /**
   * @brief varidaic initialization for multiple json files.
   * @details see ArgumentExtender::_init() for more details.
   */
  template<typename ... TS> ArgumentExtender( const std::string& first, TS ... others );

  virtual ~ArgumentExtender ();

  ArgumentExtender& AddOptions( const po::options_description& optdesc );
  void              ParseOptions( int argc, char** argv );
  bool              CheckArg( const std::string& opt ) const;

  /**
   * @brief template function for getting the user input of an option.
   * @details the user is responsible for casting the input to an appropriate
   * type.
   */
  template<typename TYPE = std::string>
  TYPE Arg( const std::string& opt ) const;

  /**
   * @brief Returning user input if exists, returning default value other wise
   * @details The user is responsible for casting the input to an appropriate
   * type.
   */
  template<typename TYPE = std::string>
  TYPE ArgOpt( const std::string& opt, const TYPE& val ) const;

  /**
   * @brief template function for getting the user input list to an option.
   * @details Notice that this must be used together with a
   * boost::program_options' multitoken when defining the options description.
   * Again the user is responsible for providing the appropriate casting.
   */
  template<typename TYPE>
  std::vector<TYPE> ArgList( const std::string& opt ) const;

  /**
   * @brief template function for getting the extend object to an options.
   * @details The user is responsible for providing the correct string tag, and
   *          providing the correct type to case the data in the json file.
   */
  template<typename TYPE>
  TYPE ArgExt( const std::string& opt, const std::string& exttag ) const;

  /** @brief Constant access to internal property tree instance. */
  inline const pt::ptree&
  Tree() const { return _exttree; }

  /**
   * @brief Constant access to internal options description instance.
   */
  inline const po::options_description&
  Description() const { return _optdesc; }

  /**
   * @brief Constant access to internal argument value map instance.
   */
  inline const po::variables_map&
  Args() const { return _argmap; }

  /**
   * @brief Simple wrapper for Argument naming scheme for filename generation.
   */
  struct ArgPathScheme
  {
    ArgPathScheme( const std::string& opt, const std::string fstring ) :
      option( opt ),
      pathstring( fstring ){}
    ArgPathScheme( const std::string& opt ) :
      option( opt ),
      pathstring( opt ){}
    std::string option;
    std::string pathstring;
  };

  typedef std::vector<ArgPathScheme> PathScheme;

  void SetFilePrefix( const fs::path );
  void SetDirScheme( const PathScheme& );
  void AddDirScheme( const ArgPathScheme& );
  void AddDirScheme( const PathScheme& );
  void SetNameScheme( const PathScheme& );
  void AddNameScheme( const ArgPathScheme& );
  void AddNameScheme( const PathScheme& );

  /**
   * @brief
   * @{
   */
  fs::path MakeFile( const std::string&, const std::string& ) const;
  fs::path MakePDFFile( const std::string& ) const;
  fs::path MakePNGFile( const std::string& ) const;
  fs::path MakeTXTFile( const std::string& ) const;
  fs::path MakeTEXFile( const std::string& ) const;
  /** @} */

  void PrintHelpAndExit() const;

protected:

  /** @brief Mutable access to internal property tree instance. */
  inline pt::ptree&
  Tree(){ return _exttree; }

  /** @brief Mutable access to internal options description instance. */
  inline po::options_description&
  Description(){ return _optdesc; }

  /** @brief Mutable access to internal argument value map instance. */
  inline po::variables_map&
  Args(){ return _argmap; }


private:
  pt::ptree _exttree;
  po::options_description _optdesc;
  po::variables_map _argmap;

  fs::path _prefix;
  PathScheme _dirscheme;
  PathScheme _namescheme;

  void        _init( const std::vector<std::string>& filelist );
  std::string genPathString( const ArgPathScheme& ) const;

  // Helper functions for path generation.
  bool IsBooleanFlag( const std::string& options ) const;
  bool IsMultiToken( const std::string& options ) const ;

  std::string genPathString_Boolean( const ArgPathScheme& ) const ;
  std::string genPathString_List( const ArgPathScheme& ) const ;
  std::string genPathString_Single( const ArgPathScheme& ) const ;
};



/*-----------------------------------------------------------------------------
 *  Template implementation
   --------------------------------------------------------------------------*/
template<typename ... TS>
ArgumentExtender::ArgumentExtender( const std::string& first, TS ... others )
{ _init( MakeVector<std::string>( first, others ... ) ); }

template<typename T>
T
ArgumentExtender::Arg( const std::string& opt ) const
{
  if( !CheckArg( opt ) ){
    throw std::invalid_argument(
      usr::fstr( "Option [%s] was not provided as a program options, "
        ", see --help output", opt ) );
  }
  return _argmap[opt].as<T>();
}

template<typename T>
T
ArgumentExtender::ArgOpt( const std::string& opt, const T& val ) const
{ return CheckArg( opt ) ? Arg<T>( opt ) : val; }

template<typename T>
std::vector<T>
ArgumentExtender::ArgList( const std::string& opt ) const
{
  if( !CheckArg( opt ) ){
    throw std::invalid_argument(
      usr::fstr( "Option [%s] was not provided as a program options, "
        ", see --help output", opt ) );
  }
  return _argmap[opt].as<std::vector<T> >();
}

template<typename T>
T
ArgumentExtender::ArgExt( const std::string& opt, const std::string& exttag ) const
{
  return usr::GetSingle<T>( Tree(), opt, Arg<std::string>( opt ), exttag );
}

}/* usr */


namespace boost
{
namespace program_options
{
/**
 * @brief Template short hand for assigning a option with a default value
 */
template<typename T>
value_semantic*
defvalue( const T def )
{
  return boost::program_options::value<T>()->default_value( def );
}

/**
 * @brief Template short had for assign a option to take multiple tokens.
 */
template<typename T>
value_semantic*
multivalue()
{
  return boost::program_options::value<std::vector<T> >()->multitoken();
}

}/* program_options */
}/* boost */

#endif/* end of include guard: USERUTILS_COMMON_OPTSNAMER_HPP */
