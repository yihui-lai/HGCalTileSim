/**
 * @file ObjectExtendedVars.hpp
 * @author [Yi-Mu "Enoch" Chen](https://github.com/yimuchen/)
 * @brief additional calculation and caching of physics variables
 */

#ifndef USERUTILS_PHYSUTILS_OBJECTEXTENDEDVARS_HPP
#define USERUTILS_PHYSUTILS_OBJECTEXTENDEDVARS_HPP

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TLorentzVector.h"

namespace usr {

/**
 * @defgroup extendedvar Extended variables.
 * @ingroup  PhysUtils
 * @brief    Handy functions store calculating commonly analysis variables
 *           not directly implemented in CMSSW
 * @{
 */

namespace mu { // muon functions, implemented in file MuonVar
extern double PfIso( const pat::Muon& );
extern double TrackIso( const pat::Muon& );
} /* mu */

/**
 * @brief Adding lorentz vector into pat::objects
 *
 * Useful when analysis requires re-evaluation of kinematics. We will be using
 * userFloats to store the numbers (rather than the userObj as that is more
 * finicky to use). If the variables used to store the momentum components
 * already exists as a userFloat, it will be overwritten without warning.
 */
template<typename T>
void
AddLorentzVector( T& obj, const TLorentzVector& vec, const std::string& name )
{
  static const bool overwrite = true;
  obj.addUserFloat( name+"pt",  vec.Pt(),     overwrite );
  obj.addUserFloat( name+"eta", vec.Eta(),    overwrite );
  obj.addUserFloat( name+"phi", vec.Phi(),    overwrite );
  obj.addUserFloat( name+"e",   vec.Energy(), overwrite );
}

/**
 * @brief Checking if pat::object has all user floats required to construct
 *        custom variables.
 */
template<typename T>
bool
HasLorentzVector( const T& obj, const std::string& name )
{
  if( obj.hasUserFloat( name+"pt" )
      && obj.hasUserFloat( name+"eta" )
      && obj.hasUserFloat( name+"phi" )
      && obj.hasUserFloat( name+"e" ) ){
    return true;
  } else {
    return false;
  }
}

/**
 * @brief Getting the internally stored lorentz vector.
 *
 * In case the variables doesn't exists for constructing a lorentz vector, the
 * four momentum of the object will be returned instead.
 */
template<class T>
TLorentzVector
GetLorentzVector( const T& obj, const std::string& name )
{
  TLorentzVector ans;
  double pt;
  double eta;
  double phi;
  double e;
  if( HasLorentzVector( obj, name ) ){
    pt  = obj.userFloat( name+"pt" );
    eta = obj.userFloat( name+"eta" );
    phi = obj.userFloat( name+"phi" );
    e   = obj.userFloat( name+"e" );
  } else {
    pt  = obj.pt();
    eta = obj.eta();
    phi = obj.phi();
    e   = obj.energy();
  }
  ans.SetPtEtaPhiE( pt, eta, phi, e );
  return ans;
}

/** @} */

}/* usr */

#endif/* end of include guard: USERUTILS_PHYSUTILS_OBJECTEXTENDEDVARS_HPP */
