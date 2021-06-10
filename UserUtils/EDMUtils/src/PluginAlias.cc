/**
 * @file    PluginAlias.cc
 * @brief   Implementation of non plugin alias functions
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 */

#include "UserUtils/EDMUtils/interface/PluginAlias.hpp"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TFile.h"
#include <string>

namespace usr {

/**
 * @brief returing the string to a path in the parameterset if an
 *        EDM::FileInPath was used.
 */
std::string
PluginAlias::GetFilePath(
 const edm::ParameterSet& config,
 const std::string& filetag )
{
  return config.getParameter<edm::FileInPath>( filetag ).fullPath();
}

/**
 * @brief Given a file path in the form of an EDM::FileInPath, and a object tag
 *        as a parameter set string. This function return a clone to the
 *        object stored in a file for analysis use.
 */
TObject*
PluginAlias::GetFileObj(
  const edm::ParameterSet& config,
  const std::string&       filetag,
  const std::string&       objtag )
{
  const std::string filename = GetFilePath( config, filetag );
  const std::string objname  = config.getParameter<std::string>( objtag );
  TFile* file                = TFile::Open( filename.c_str() );
  TObject* ans               = file->Get( objname.c_str() )->Clone();
  file->Close();
  return ans;
}

}/* usr */
