/**
 * @file    TriggerMatching.hpp
 * @brief   Trigger object matching functions.
 * @author  [Yi-Mu "Enoch" Chen](https://github.com/yimuchen)
 *
 *  @details Attempts to match a reco::Candidate object to a trigger object in
 *  a list of trigger objects with a specified trigger within a given cone
 *  size (defaults to 0.1).
 *
 *  The triggerNames object required should be  extracted directly from
 *  the edm::Event and edm::TriggerResults class ese .
 *  The requested trigger (reqtrigger) and requested label (reqlabel) could
 *  have globbing syntax (e.x. HLT_IsoMu27_v* for multiple versions)
 */

#ifndef USERUTILS_PHYS_UTILS_HPP
#define USERUTILS_PHYS_UTILS_HPP

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <string>

#include "TLorentzVector.h"

namespace usr {


/**
 * @defgroup triggermatching Trigger Matching
 * @brief Functions for performing trigger matching.
 * @ingroup PhysUtils
 *
 * @{
 */
extern bool HasTriggerMatch(
  const reco::Candidate&                        recoobj,
  const pat::TriggerObjectStandAloneCollection& trgobjlist,
  const std::string&                            reqtrigger,
  const std::string&                            reqlabel,
  const edm::TriggerNames&                      trgnames,
  const double                                  mindelta_r = 0.1
  );

extern TLorentzVector TriggerMatchP4(
  const reco::Candidate&                        recoobj,
  const pat::TriggerObjectStandAloneCollection& trgobjlist,
  const std::string&                            reqtrigger,
  const std::string&                            reqlabel,
  const edm::TriggerNames&                      trgnames,
  const double                                  mindelta_r = 0.1
  );


/*-----------------------------------------------------------------------------
 *  MatchTrgObj
 *   Return a vector::const_iterator to the first tracker object matched with
 *   reco::Candidate object within a radius of midelta_r. returns the end()
 *   iterator if match is not found.
 *   It is not recommended to use this function directly, but to use by the
 *   functions above.
   --------------------------------------------------------------------------*/
extern pat::TriggerObjectStandAloneCollection::const_iterator MatchTrgObj(
  const reco::Candidate&                        recoobj,
  const pat::TriggerObjectStandAloneCollection& trgobjlist,
  const std::string&                            reqtrigger,
  const std::string&                            reqfilter,
  const edm::TriggerNames&                      trgnames,
  const double                                  mindelta_r = 0.1
  );

/*-----------------------------------------------------------------------------
 *  TrigObjMatchPathFilter
 *   Check if single trigger object maches desired path and filter label.
 *   Input could be glob notation and will return first glob match results
   --------------------------------------------------------------------------*/
extern bool TrigObjMatchPathFilter(
  const pat::TriggerObjectStandAlone& obj,
  const std::string&                  path,
  const std::string&                  filter
  );

/** @} */

}/* usr */

#endif/* end of include guard: USERUTILS_PHYS_UTILS_HPP */
