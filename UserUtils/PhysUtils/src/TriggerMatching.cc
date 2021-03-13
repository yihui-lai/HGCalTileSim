/**
 * @file TriggerMatching.cc
 * @brief Implementing trigger matching recipes
 * @author [Yi-Mu "Enoch" Chen](https://github.com/yimuchen/)
 */
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "UserUtils/PhysUtils/interface/TriggerMatching.hpp"
#include "UserUtils/Common/interface/STLUtils/StringUtils.hpp"

#include <regex>

#include "TLorentzVector.h"

/**
 * @defgroup triggermatching Trigger Matching
 * @details
 * @details
 * Main reference could be found in the documentation of
 * [general trigger objects use](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger).
 *
 * In general, you will need provide :
 * - physical object you want to perform the matching (ex, a pat::Jet,
 *   or pat::Muon).
 * - The list of trigger objects stored in the EDM file (should be under the
 *   input tag of "selectedPatTrigger" or something similar).
 * - The trigger bits in the EDM file (should be under the input tag
 *   "TriggerResults::HLT" or similar).
 * - Trigger paths and filter label. For paths, please consult your physics
 *   analysis group for recommendations. For the filter label of objects within
 *   a trigger, please consult the [trigger database](https://cmsweb.cern.ch/confdb/)
 *   (warning! This website is notoriously slow!), or get the trigger
 *   configuration via the [`hltGetConfiguration` command](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT).
 */

namespace usr {

/**
 * @brief Check if a reco::Candidate object has any matching trigger objects
 *        in the provided trigger list with the specific trigger name and
 *        filter withing the specific \f$\Delta R\f$ requirements.
 *
 * This function is simple used for checking the exists, with no intermediates
 * results stored. If you need to store the physical properties of the matching
 * trigger objects, consider using the other functions.
 */
bool
HasTriggerMatch(
  const reco::Candidate&                        recoobj,
  const pat::TriggerObjectStandAloneCollection& trgobjlist,
  const std::string&                            reqtrigger,
  const std::string&                            reqfilter,
  const edm::TriggerNames&                      trgnames,
  const double                                  mindelta_r
  )
{
  auto match = MatchTrgObj(
    recoobj,
    trgobjlist,
    reqtrigger,
    reqfilter,
    trgnames,
    mindelta_r );
  return match != trgobjlist.end();
}

/**
 * @brief Returning the four momentum of the first matched trigger object.
 *
 * In case a trigger object is not found, the function will return a
 * (0,0,0,0) four momentum.
 */
TLorentzVector
TriggerMatchP4(
  const reco::Candidate&                        recoobj,
  const pat::TriggerObjectStandAloneCollection& trgobjlist,
  const std::string&                            reqtrigger,
  const std::string&                            reqfilter,
  const edm::TriggerNames&                      trgnames,
  const double                                  mindelta_r
  )
{
  auto match = MatchTrgObj(
    recoobj,
    trgobjlist,
    reqtrigger,
    reqfilter,
    trgnames,
    mindelta_r );
  if( match != trgobjlist.end() ){
    return TLorentzVector( match->px(), match->py(), match->pz(), match->energy() );
  } else {
    return TLorentzVector( 0, 0, 0, 0 );
  }
}

/**
 * @brief Check if a reco::Candidate object has any matching trigger objects
 *        in the provided trigger list with the specific trigger name and
 *        filter withing the specific \f$\Delta R\f$ requirements.
 *
 * This function returns the iterator to the trigger object, so the user
 * is responsible for checking if it is out of bounds of the original trigger
 * object list. The first object would be returned, so if you are checking
 * for the existence of multiple objects, you might need write your own
 * function.
 */
pat::TriggerObjectStandAloneCollection::const_iterator
MatchTrgObj(
  const reco::Candidate&                        recoobj,
  const pat::TriggerObjectStandAloneCollection& trgobjlist,
  const std::string&                            reqtrigger,
  const std::string&                            reqfilter,
  const edm::TriggerNames&                      trgnames,
  const double                                  mindelta_r
  )
{
  for( auto objiter = trgobjlist.begin(); objiter != trgobjlist.end(); ++objiter ){
    auto obj = *objiter;  // making duplicate
    obj.unpackPathNames( trgnames );

    if( TrigObjMatchPathFilter( obj, reqtrigger, reqfilter )
        && reco::deltaR( recoobj, obj ) < mindelta_r ){
      return objiter;
    }
  }

  return trgobjlist.end();
}

/**
 * @brief Given a trigger object, check if it has trigger paths and
 *       trigger filters matching the input.
 *
 * Globbing notations are allowed for both the trigger path and filter label
 * input (i.e. given the path "HLT_Ele27_v*", the function will consider
 * the path matched if any version of the paths used for the trigger object.)
 * Only if both the path and filter label are found to be matched will the
 * function return true.
 */
bool
TrigObjMatchPathFilter(
  const pat::TriggerObjectStandAlone& obj,
  const std::string&                  path,
  const std::string&                  filter
  )
{
  bool hasmatchpath   = false;
  bool hasmatchfilter = false;

  const std::regex pathreg( usr::GlobToRegex( path ) );

  for( const auto& objpath : obj.pathNames() ){
    if( std::regex_match( objpath, pathreg )
        && obj.hasPathName( objpath ) ){
      hasmatchpath = true;
      break;
    }
  }

  const std::regex filterreg( usr::GlobToRegex( filter ) );

  for( const auto& objfilter : obj.filterLabels() ){
    if( std::regex_match( objfilter, filterreg )
        && obj.hasFilterLabel( objfilter ) ){
      hasmatchfilter = true;
      break;
    }
  }

  return hasmatchpath && hasmatchfilter;
}

} /* usr */
