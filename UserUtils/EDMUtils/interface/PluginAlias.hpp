/**
 * @file    PluginAlias.hpp
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 * @brief   Defining alias functions for EDM plugins via intermediate virtual
 *          class
 */

#ifndef USERUTILS_EDMUTILS_PLUGINALIAS_HPP
#define USERUTILS_EDMUTILS_PLUGINALIAS_HPP

#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TObject.h"

#include <string>

namespace usr {

/**
 * @brief class for reducing the verbosity of the various edm plugin calls.
 * @ingroup EDMUtils
 *
 * @details One I find immensely annoying to code is the `consume` calls to
 *          create an EDMToken object. What makes it worse is that the consumes
 *          call is an protected member function of the edm plugin's base class.
 *          Writing an alias to the function, therefore, requires all future
 *          plugins to inherit from this alias class.
 */
class PluginAlias : public virtual edm::EDConsumerBase
{
public:

  /**
   * Like all EDM plugins, this requires a parameterset to initialize. A
   * reference to the parameterset is stored for less verbose parsing.
   */
  PluginAlias( const edm::ParameterSet& config ) : _config( config ){}
  virtual
  ~PluginAlias(){}

protected:

  /**
   * @brief template class for getting the EdmToken from the initializer
   *        configuration parameter set.
   *
   * @tparam TYPE  the type of object attempting to be consumed must be user
   *               specified.
   */
  template<typename TYPE>
  inline edm::EDGetToken
  GetToken( const std::string& configtag )
  {
    return consumes<TYPE>( _config.getParameter<edm::InputTag>( configtag ) );
  }

  /**
   * @brief Same template function as GetToken() except for getting Run level
   *        objects.
   */
  template<typename TYPE>
  inline edm::EDGetToken
  GetRunToken( const std::string& configtag )
  {
    return consumes<TYPE, edm::InRun>(
      _config.getParameter<edm::InputTag>( configtag ) );
  }

  /**
   * @brief getting a clone of an object described in a file. See static
   *        function for full documentation.
   */
  inline TObject*
  GetFileObj( const std::string& filetag, const std::string objtag )
  { return GetFileObj( _config, filetag, objtag ); }

  /**
   * @brief getting the full path to a file. See static function for full
   *        documentation.
   */
  inline std::string
  GetFilePath( const std::string& filetag )
  { return GetFilePath( _config, filetag ); }

  static std::string GetFilePath(
    const edm::ParameterSet&,
    const std::string& );

  static TObject* GetFileObj(
    const edm::ParameterSet& config,
    const std::string&       file,
    const std::string&       objname );

private:
  const edm::ParameterSet& _config;

};

}/* usr */

#endif/* end of include guard: USERUTILS_EDMUTILS_PLUGINALIAS_HPP */
