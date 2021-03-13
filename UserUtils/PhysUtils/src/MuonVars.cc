/**
 * @file MuonVars.cc
 * @brief Muon variables calculation.
 * @author [Yi-Mu "Enoch" Chen](https://github.com/yimuchen/)
 */
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <algorithm>

namespace usr {

namespace mu {

/**
 * @addtogroup extendedvar
 * @{
 */

/**
 * @brief The relative particle-flow isolation
 *
 * Full documentation of the variable [here](https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Trigger).
 * Since the cut values are every changing, this will not be provided here.
 */
double
PfIso( const pat::Muon& mu )
{
  const double chghad = mu.pfIsolationR04().sumChargedHadronPt;
  const double neuhad = mu.pfIsolationR04().sumNeutralHadronEt;
  const double pho    = mu.pfIsolationR04().sumPhotonEt;
  const double pu     = mu.pfIsolationR04().sumPUPt;
  const double pt     = mu.pt();

  return ( chghad+std::max( 0.0, neuhad+pho-0.5*pu ) )/pt;
}

/**
 * @brief Relative track isolation
 *
 * Full documentation of the variable [here](https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Trigger).
 * Since the cut values are every changing, this will not be provided here.
 */
double
TrackIso( const pat::Muon& mu )
{
  return mu.isolationR03().sumPt / mu.pt();
}

/** @} */

} /* mu */

} /* usr */
